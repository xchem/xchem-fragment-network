# NodeProvider
#
# SocketIO:  Note that NodeProvider connects to the ChemiReg server using a WebSocket and SocketIO.
# The NodeProvider class below contains only the parts of the ChemiReg server protocol required for AuditClass.
# In the future a transpiled version of the SATurn class saturn.db.NodeProvider will be used instead but for that
# to happen saturn.client.core.ClientCore needs to be made Python compatible.

from SATurn import haxe_Serializer
from SATurn import haxe_ds_StringMap
from SATurn import _hx_AnonObject

import urllib.request
import urllib.parse
import json
import ijson

import logging
import time

import os
import base64
from astropy import log

from socketIO_client_nexus import SocketIO, LoggingNamespace

logging.getLogger('socketIO-client').setLevel(logging.INFO)
logging.basicConfig()

class NodeProvider(object):
    def __init__(self, hostname, port, username, password, cb):
        self.hostname = hostname
        self.port = port
        self.username = username
        self.password = password
        self.auth_token = None
        self.socketIO = None

        self.next_msg_id = 0

        self.msg_id_to_callback = {}

        self._after_connect = cb

    def login(self):
        login_url = self.hostname + ':' + str(self.port) + '/login'

        values = {
            'username': self.username, 'password': self.password
        }

        data = urllib.parse.urlencode(values).encode('ascii')

        request = urllib.request.Request(login_url, data)

        print('Connecting to ' + login_url)

        while True:
            try:
                with urllib.request.urlopen(request) as response:
                    content = response.read()

                    content_obj = json.loads(content.decode('ascii'))

                    if 'token' in content_obj:
                        self.auth_token = content_obj['token']

                        self.configure_socket()
                        break
                    else:
                        raise Exception('Authentication failed no token in response')
            except urllib.error.URLError as e:
                print('Sleeping!')
                print(e)
                time.sleep(10)

    def configure_socket(self):

        self.socketIO = SocketIO('https://globalchemireg.sgc.ox.ac.uk', self.port, LoggingNamespace,
                                 transports=['websocket'])
        time.sleep(3)
        self.socketIO.on('authenticated', self._socket_authenticated)
        self.socketIO.emit('authenticate', {'token': self.auth_token})
        self.socketIO.on('__response__', self._process_response)

        self.socketIO.on('reconnect', self._socket_authenticated)
        self.socketIO.on('connect', self._socket_authenticated)
        self.socketIO.on('disconnect', self._disconnected)
        self.socketIO.on('open', self._socket_authenticated)

        self.socketIO.wait()

    def _disconnected(self):
        print('Disconnected')

    def _process_response(self, data):
        if 'bioinfJobId' in data:
            msg_id = data['bioinfJobId']

            if msg_id in self.msg_id_to_callback:
                error = None
                json = None

                if 'error' in data:
                    error = data['error']
                elif 'json' in data and 'error' in data['json']:
                    error = data['json']['error']
                else:
                    error = None

                if 'json' in data:
                    json = data['json']

                self.msg_id_to_callback[msg_id](error, json)

                self.msg_id_to_callback.pop(msg_id)
            else:
                print('Warning message not found ' + msg_id)

    def handle_fetch(self, error, data):
        print(data)

    def _authenticate_socket(self):
        print('Authenticating')
        self.socketIO.emit('authenticate', {token: self.auth_token})

    def _socket_authenticated(self):
        self.after_connect()

    def after_connect(self):
        self._after_connect()

    def set_after_connect(self, cb):
        self._after_connect = cb

    def get_by_named_query(self, query_id, data, cb):
        json = {"queryId": query_id, 'parameters': self._serialise(data)}

        def _cb(error, json):
            cb(error, data['objects'])

        self.run_query('_remote_provider_._data_request_objects_namedquery', json, cb)

    def upload_file(self, filename, cb):
        f = open(filename, 'rb')

        chunk_size = 1024 * 60000

        scope_obj = {'file_identifier': None}

        def upload_chunk():
            byte_buf = f.read(chunk_size)

            eof = len(byte_buf) != chunk_size

            def _cb(error, upload_id):
                if error is not None:
                    f.close()

                    cb(error, None)

                scope_obj['file_identifier'] = upload_id

                if not eof:
                    upload_chunk()
                else:
                    f.close()

                    cb(None, scope_obj['file_identifier'])

            self.upload_bytes_as_file(byte_buf, scope_obj['file_identifier'], _cb)

        upload_chunk()

    def upload_bytes_as_file(self, contents, file_identifer, cb):
        contents_b64 = base64.b64encode(contents).decode('ascii')
        json = {'contents': contents_b64, 'file_identifier': file_identifer}

        def _cb(error, json):
            cb(error, json['upload_id'])

        self.run_query('_remote_provider_._data_request_upload_file', json, _cb)

    def _convert_dictionaries(self, obj):
        for key in obj.keys():
            if type(obj[key]) == dict:
                obj[key] = self._convert_dictionaries(obj[key])

        return _hx_AnonObject(obj)

    def _serialise(self, params):
        a = self._convert_dictionaries(params)

        param_str = haxe_Serializer.run([a])

        return param_str

    def run_query(self, api_command, json, cb):
        msg_id = self.increment_next_id()

        json['msgId'] = msg_id

        self.register_callback(msg_id, cb)

        self.socketIO.emit(api_command, json)

    def register_callback(self, msg_id, cb):
        self.msg_id_to_callback[msg_id] = cb

    def increment_next_id(self):
        i = self.next_msg_id

        self.next_msg_id += 1

        return str(i)