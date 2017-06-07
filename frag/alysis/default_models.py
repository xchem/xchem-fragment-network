from frag.alysis.template import TemplateClass
from frag.utils import parse_waters

class WaterCluster(TemplateClass):
    """
    A class to cluster waters
    """

    def __init__(self, input_raw_data):
        self.input_raw_data = input_raw_data

    def parse_data(self):
        parse_waters(self.input_raw_data)
