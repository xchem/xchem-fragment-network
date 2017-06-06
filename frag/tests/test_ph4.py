from frag.utils import RDKitPh4
import unittest
from rdkit import Chem


class Ph4Test(unittest.TestCase):
    sdf_data = """
     RDKit          3D

  8  8  0  0  0  0  0  0  0  0999 V2000
   -0.2375    1.2479   -0.2221 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5033    1.0878    0.3249 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9231   -0.1960    0.6598 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4265   -1.1855   -0.1809 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0741   -1.1112   -0.4792 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5447    0.1297   -0.4205 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0371    0.2224   -0.4859 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5826   -0.1950    0.8038 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  1  0
  6  1  1  0
M  END

"""

    def test_ph4_parser(self):
        """
        Test the pharmacophore generation
        :return:
        """
        rdmol = Chem.MolFromMolBlock(self.sdf_data)
        rdkit_ph4 = RDKitPh4()
        feats = rdkit_ph4.generate_ph4_for_mol(rdmol=rdmol)
        self.assertEqual(len(feats),5)


