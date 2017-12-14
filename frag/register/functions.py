from .chemireg import ChemiReg
import pprint

HOST_NAME = 'https://globalchemireg.sgc.ox.ac.uk'
PORT = 443

class ChemiRegInterface():

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.client.close()

    def __enter__(self):
        self.client.connect()
        return self

    def __init__(self,user_name, password, project_name="OxXChem"):
        self.client = ChemiReg(HOST_NAME, PORT, user_name, password)
        self.pp = pprint.PrettyPrinter(indent=4)
        self.project_name=project_name

    def _handle_query(self,query):
        self.pp.pprint(query)
        if query['error'] is not None:
            raise Exception(query['error'])
        return query['objects']

    def _handle_update_input(self, data_to_update):
        """
        Handle the update input so that it contains ids
        :param data_to_update: the input data needing updating
        :return: the dictionary with ids as the key and the data as the value
        """
        data_with_ids = [x for x in data_to_update if 'id' in x]
        if len(data_with_ids) == 0:
            print("Warning - no entries to update")
            return False
        data_not_updated = len(data_to_update) - len(data_with_ids)
        if data_not_updated != 0:
            print("Warning - "+str(data_not_updated)+" datapoints cannot be update.")
        id_as_key = {}
        for x in data_with_ids:
            id_as_key[x['id']] = x
        return id_as_key,list(id_as_key.keys())

    def _update_data_on_ids(self,query,ids_as_key):
        """
        Update the data on a query based on data in a corresponding dictoonary
        :param query:
        :param ids_as_key:
        :return:
        """
        changes = []
        for compound in query['objects']:
            data = ids_as_key[compound['id']]
            for key in data:
                compound[key] = data[key]
            changes.append(compound)
        return changes

    def update_project_name(self,project_name):
        self.project_name=project_name


    def fetch_ids(self, ids_to_fetch):
        """
        Fetch ids
        :param ids_to_fetch:
        :param project_name:
        :return:
        """
        query = self.client.fetch(ids_to_fetch, self.project_name, 0, len(ids_to_fetch))
        return self._handle_query(query)

    def update_ids(self,data_to_update):
        """

        :param data_to_update: list of dictionaries
                       data = [
            {
                'id': 'DI000016a',
                'supplier_id': 'ENAMINE',
            },
            {
                'id': 'DI000016a',
                'supplier_id': 'ENAMINE',
            }
        :return:
        """
        ids_as_key,ids_to_update = self._handle_update_input(data_to_update)
        if not ids_as_key:
            return
        # Now perform the query
        query = self.client.fetch(ids_to_update, self.project_name, 0, len(ids_to_update))
        changes = self._update_data_on_ids(query,ids_as_key)
        query = self.client.save_changes(changes, self.project_name)
        return self._handle_query(query)

    def insert_json(self, input_data):
        """
        In
        :param input_data: JSON to save of the form
                data = [
            {
                'supplier': 'Test',
                'supplier_id': 'A',
                'classification': 'DI',
                'smiles': 'c1ccccc1'
            },
            {
                'supplier': 'Test',
                'supplier_id': 'B',
                'classification': 'DI',
                'smiles': 'c1ccccc1'
            }
        ]
        :return: the query to be scrutinised
        """
        query = self.client.save_changes(input_data, self.project_name)
        return self._handle_query(query)
