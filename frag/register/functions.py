import requests


def put_compounds(user_name, password, input_data):
    """
    Put a series of compounds up to the Chemireg API
    :param user_name: your chemireg username
    :param password: your chemireg password
    :param input_data: the JSON to upload
    :return: None
    """
    url = "https://globalchemireg.sgc.ox.ac.uk/static/chemireg/upload/"
    request = requests.post(url, json=input_data, auth=(user_name, password))
    if request.status_code == requests.codes.ok:
        return request.json()
    else:
        raise IOError
