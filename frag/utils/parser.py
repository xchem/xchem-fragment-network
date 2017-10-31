# Series of functions to parse input files
from rdkit import Chem
from frag.alysis.models import Cluster_Things,Object,Owner

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


def _parse_ligand_sdf(input_file):
    """
    Function to parse a series of ligands - return RDKit mols.
    :param input_file: the file to parse
    :param input_type: the type of ligands
    :return: the molecules parsed
    """
    mols = Chem.SDMolSupplier(input_file)
    return mols

def parse_ligands(input_file,input_type="sdf"):
    mols = _parse_ligand_sdf(input_file=input_file)
    # Now return them with their name and centre of mass


def parse_ligand_ph4s(input_file, input_type="sdf"):
    """
    Function to parse a series of ligands - return RDKit mols.
    :param input_file: the file to parse
    :param input_type: the type of ligands
    :return: the molecules parsed
    """
    mols = Chem.SDMolSupplier(input_file)
    return mols

def parse_waters(input_pdbs, input_mol=None, max_dist=10.0):
    """
    Function to parse a series of waters - return waters.
    :param input_pdb: the input PDB files
    :param input_mol: the input molecule (to use as a reference around which to limit)
    :return: tuple threes of coordinates of the waters
    """
    owner_list = []
    # First just get the waters from the file
    for input_pdb in input_pdbs:
        waters = _get_waters(open(input_pdb).readlines())
        water_coords = _get_water_coords(waters)
        out_l = []
        for water in water_coords:
            water = Object(water,"water")
            out_l.append(water)
        owner = Owner(out_l,input_pdb)
        owner_list.append(owner)
    return owner_list


### TODO Residues
def _get_res(data):
    """Helper function to get the coordinates from a load of waters."""
    rd_waters = Chem.MolFromPDBBlock("\n".join(data))
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
            out_list.append((cp.x, cp.y, cp.z))
    return out_list


def parse_residues(input_pdbs, input_mol=None, max_dist=10.0):
    """
    Function to parse a series of PDB files of proteins.
    :param input_pdb: the input PDB files - with identifiers
    :return: a dict (Key Residue -> value list of molecules)
    """
    owner_list = []
    res_dict = {}
    # First just get the waters from the file
    for input_pdb in input_pdbs:
        # Get a list of residues with RDKit mols (RDMol,RES_NAME)
        res_list = _get_res(open(input_pdb).readlines())
        for res in res_list:
            res_name = res[1]
            rd_mol = res[2]
            if res in res_dict:
                res_dict[res].append(res)
            else:
                res_dict[res] = [res]
    for res in res_dict:
        rmsd_coords = _get_res_rmsds(res_dict[res])
        out_l = []
        res = Object(rmsd_coords,res)
        out_l.append(res)
        owner = Owner(out_l,input_pdb)
        owner_list.append(owner)
    return owner_list

def get_file(file_path,output_format,file_counter):
    if output_format=="smi":
        return Chem.SmilesWriter(file_path+"_"+str(file_counter)+".smi")
    else:
        return Chem.SDWriter(file_path+"_"+str(file_counter)+".sdf")