# Series of functions to parse input files
from rdkit import Chem


def _get_waters(lines):
    """Helper function to extract waters from a PDB file"""
    return [line for line in lines]# if line[17:20] == "HOH"]


def _get_water_coords(waters):
    """Helper function to get the coordinates from a load of waters."""
    rd_waters = Chem.MolFromPDBBlock("\n".join(waters))
    out_list = []
    if rd_waters is None:
        print("Warning - unable to parse waters.")
    if rd_waters is not None:
        # Check the waters exist
        conf = rd_waters.GetConformer()
        # Delete them for this protein
        for i in range(rd_waters.GetNumAtoms()):
            cp = conf.GetAtomPosition(i)
            if rd_waters.GetAtomWithIdx(i).GetSmarts() != "O":
                print("Warning - skipping a water")
                continue
            out_list.append((cp.x,cp.y,cp.z))
    return out_list


def parse_ligands(input_file, input_type="sdf"):
    """
    Function to parse a series of ligands - return RDKit mols.
    :param input_file: the file to parse
    :param input_type: the type of ligands
    :return: the molecules parsed
    """
    mols = Chem.SDMolSupplier(input_file)
    return mols


def parse_waters(input_pdb, input_mol=None, max_dist=10.0):
    """
    Function to parse a series of waters - return waters.
    :param input_pdb: the input PDB files
    :param input_mol: the input molecule (to use as a reference)
    :return: tuple threes of coordinates of the waters
    """
    # First just get the waters from the file
    waters = _get_waters(open(input_pdb).readlines())
    water_coords = _get_water_coords(waters)
    return  water_coords


### TODO Residues

def parse_residues(input_pdb, input_mol=None, max_dist=10.0):
    """
    Function to parse a series of PDB files of proteins.
    :param input_pdb: the input PDB files
    :return: a dict (Key Residue -> value list of molecules)
    """
    # First just get the waters from the file
    waters = _get_waters(open(input_pdb).readlines())
    water_coords = _get_water_coords(waters)
    return  water_coords