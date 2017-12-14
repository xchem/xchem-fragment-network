

class Hypothesis(object):
    """
    Class for storing a basic hypothesies
    """

    def __init__(self, target_protein, base_mol_ids, driver,
                 inform_mol_ids=None, comment=None, structural_data=None):
        """
        Instantiate a new register instance for a given compound
        :param target_protein: the target protein we're working with
        :param base_mol_ids: the list of molecule ids related to this
        :param driver: the driver for this register (chemical, structural)
        :param inform_mol_ids: other molecules that inform this hpyothesis
        :param comment: an optional comment on the register
        :param structural_data: an optional structural display for an observation
        """
        self.target_protein = target_protein
        self.base_mol_ids = base_mol_ids
        self.driver = driver
        self.inform_mol_ids = inform_mol_ids
        self.comment = comment
        self.structural_data = structural_data









