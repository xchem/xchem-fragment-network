

class TemplateObject(object):

    input_type = None
    input_features = None

    def __init__(self, input_type, input_features):
        self.input_features = input_features
        self.input_type = input_type


class TemplateCluster(object):

    output_type = None
    output_centre = None
    def __init__(self, output_type, output_centre):
        self.output_type = output_type
        self.output_centre = output_centre


class TemplateClass(object):

    # The raw data source
    input_raw_data = None
    # A list of TemplateObject objects
    input_data = None
    # A list of TemplateClusters
    output_data = None

    def parse_data(self):
        """
        Function to parse data from your input format.
        :return: a list of TemplateObject
        """
        raise NotImplementedError

    def cluster_data(self):
        """
        Uses the input data to create clusters.
        :return: a list of cluster objects
        """
        raise NotImplementedError
