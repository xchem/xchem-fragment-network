# Generate ph4s for a molecule
from rdkit.Chem import ChemicalFeatures
import os



class RDKitPh4(object):

    factory = None

    def get_factory(self):
        """
        Generate the Ph4 feature factory
        :return:
        """
        if self.factory is None:
            this_dir, this_filename = os.path.split(__file__)
            data_path = os.path.join(this_dir, "data", "RDKitPh4.fdef")
            self.factory = ChemicalFeatures.BuildFeatureFactory(data_path)
        return self.factory


    def generate_ph4_for_mol(self, rdmol):
        """
        Generate a pharmacophore from an input molecule and a feature factory.
        :param rdmol: the input RDKit molecule
        :param factory: the feature factory
        :return: a list of 4 tuples (x,y,z, feature)
        """
        feats = self.get_factory().GetFeaturesForMol(rdmol)
        return [(feat.GetPos().x,feat.GetPos().y,feat.GetPos().z,feat.GetType()) for feat in feats]