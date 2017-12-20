from .nodeprovider import NodeProvider

from threading import Lock
import threading
import os
import pprint
import uuid
from time import sleep

class ChemiReg(object):
    def __init__(self, hostname, port, username, password):
        # ChemiReg hostname include protocol (https/http)
        self.hostname = hostname

        # ChemiReg port
        self.port = port

        # ChemiReg username
        self.username = username

        # ChemiReg password
        self.password = password

        # Lock to synchronize threads
        self.lock = Lock()

        # Variable to store web-service call results
        self.blocking_results = {'ready': False}

        # Variable to instruct web-service thread
        self.blocking_input = {}

        # Used to communicate between threads
        self.wait_event = threading.Event()

        # Thread for listen for messages
        self.web_socket_thread = None

    def set_ready(self, ready):
        self.lock.acquire()
        try:
            self.blocking_results['ready'] = ready
        finally:
            self.lock.release()

    def set_state(self, error, objects, ready):
        self.lock.acquire()
        try:
            self.blocking_results['error'] = error
            self.blocking_results['objects'] = objects
            self.blocking_results['ready'] = ready
        finally:
            self.lock.release()

    def is_ready(self):
        self.lock.acquire()
        try:
            return self.blocking_results['ready']
        finally:
            self.lock.release()  # release lock, no matter what

    def block_call(self):
        while True:
            if self.is_ready():
                return self.blocking_results
            else:

                sleep(0.05)

    def _process_response(self, error, json):
        if json['objects'] == None and json['error']:
            objs = {}
        elif 'objects' in json and 'upload_set' in json['objects'][0]:
            objs = json['objects'][0]['upload_set']
        elif 'objects' in json and 'refreshed_objects' in json['objects'][0]:
            objs = json['objects'][0]['refreshed_objects']
        else:
            objs = json['objects'][0]

        self.set_state(error, objs, True)

    def connect(self):
        # Initialise NodeProvider which handles communication with ChemiReg via a WebSocket
        self.provider = NodeProvider(self.hostname, self.port, self.username, self.password, None)

        # This function is called after the connection has been established
        def after_connect():
            # Used to inform the initiating thread that we are connected (releases the block)
            self.set_ready(True)

        def start_login():
            self.provider._after_connect = after_connect
            self.provider.login()

        self.web_socket_thread = threading.Thread(target=start_login)
        self.web_socket_thread.start()

        while True:
            if self.is_ready():
                break


    def _close(self):
        self.provider.socketIO._close()

    def close(self):
        self._close()

        self.web_socket_thread.join()

    def register_sdf(self, sdf, project_name, config, cb):
        self.set_ready(False)

        def _process_upload(error, upload_key):
            self.set_state(error, upload_key, True)
            
        self.provider.upload_file(sdf, _process_upload)

        self.block_call()

        json = {
            'upload_key_sdf': self.blocking_results['objects'],
            'name': os.path.basename(sdf),
            '_username': None,
            'project_name': project_name,
            'upload_defaults': config
        }

        self.set_ready(False)

        def _process_registration(error, obj):
            upload_key = None

            if error is None:
                upload_key = obj['objects'][0]['upload_id']

            self.set_state(error, upload_key, True)

        self.provider.get_by_named_query('saturn.db.provider.hooks.ExternalJsonHook:SDFRegister', json, _process_registration)

        return self.block_call()


    def save_changes(self, objects, project):
        changes = {}

        next_id = -1

        for object in objects:
            if 'id' in object:
                changes[object['id']] = object
            else:
                changes[str(next_id)] = object

                next_id -= 1

        self.set_ready(False)

        json = {
            'project_name': project,
            'save_changes': self.provider._convert_dictionaries(changes),
            '_username': None
        }

        self.provider.get_by_named_query('saturn.db.provider.hooks.ExternalJsonHook:SDFRegister',json, self._process_response)

        return self.block_call()

    def query(self,command_name, arguments):
        self.set_ready(False)

        self.provider.get_by_named_query(command_name, arguments, self._process_response)

        return self.block_call()

    def fetch(self, ids, project, from_row, to_row):
        arguments = {
            'action': 'fetch_exact',
            'task': 'fetch',
            'ids': ids,
            '_username': None,
            'project_name': project,
            'from_row': from_row,
            'to_row': to_row
        }

        return self.query('saturn.db.provider.hooks.ExternalJsonHook:Fetch',arguments)

    def fetch_count(self, ids, project, from_row, to_row, cb):
        arguments = {
            'action': 'fetch_exact',
            'task': 'count',
            'ids': ids,
            '_username': None,
            'project_name': project,
            'from_row': from_row,
            'to_row': to_row
        }

        return self.query('saturn.db.provider.hooks.ExternalJsonHook:Fetch', arguments)

    def fetch_upload_set(self, upload_id, project, from_row, to_row, cb):
        arguments = {
            'action': 'fetch_upload',
            'task': 'fetch',
            'upload_id': upload_id,
            '_username': None,
            'project': project,
            'from_row': from_row,
            'to_row': to_row
        }

        return self.query('saturn.db.provider.hooks.ExternalJsonHook:Fetch', arguments)

class ChemiRegTests(object):
    def __init__(self, username, password, host='https://globalchemireg.sgc.ox.ac.uk', port=443):
        self.client = ChemiReg(host, port, username, password) # type : ChemiReg
        self.client.connect()
        self.pp = pprint.PrettyPrinter(indent=4)

    def run_tests(self):
        self.fetch_test1()
        self.insert_json_test1()
        self.update_test1()

        self.after_tests()

    def fetch_test1(self):
        query = self.client.fetch(['NU000767a'], 'SGC - Oxford', 0, 1)

        self.pp.pprint(query)

        assert query['error'] is None
        assert len(query['objects']) == 1

        print('Fetch Test 1 - Passed')

    def update_test1(self):
        query = self.client.fetch(['DI000016a'], 'OxXChem', 0, 1)

        assert query['error'] is None
        assert len(query['objects']) == 1

        compound = query['objects'][0]

        changes = [{
            'id': str(compound['id']),
            'supplier_id':str(uuid.uuid4())
        }]

        query = self.client.save_changes(changes, 'OxXChem')

        self.pp.pprint(query)

        assert query['error'] is None

        saved_objects = query['objects']

        assert len(saved_objects.keys()) == 1

        assert list(saved_objects.values())[0]['supplier_id'] == changes[0]['supplier_id']

    def insert_json_test1(self):
        changes = [
            {
                'supplier': 'Test',
                'supplier_id': 'A',
                'classification': 'DI'
            },
            {
                'supplier': 'Test',
                'supplier_id': 'B',
                'classification': 'DI'
            }
        ]

        query = self.client.save_changes(changes, 'OxXChem')

        self.pp.pprint(query)

        assert query['error'] is None

        saved_objects = query['objects']

        assert len(saved_objects.keys()) == 2

    def after_tests(self):
        self.client.close()

