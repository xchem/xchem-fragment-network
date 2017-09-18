from rdkit import Chem

def get_num_ring_atoms(input_mol):
    """
    Get the number of ring atoms
    :param input_mol:
    :return:
    """
    num_ring_atoms = 0
    split_index =0
    split_indices = []
    for atom in input_mol.GetAtoms():
        if atom.IsInRing():
            num_ring_atoms+=1
        else:
            split_indices.append(split_index)
        split_index += 1
    return num_ring_atoms,split_indices


def simplified_graph(input_smiles):
    """
    The simplified graph representation for the edge uses the daylight
function dt_molgraph to
1) set all bond orders to one, 2) remove aromaticity, 3??"set the hydrogen count,
4) remove charges and set masses to zero. 5)  Additionally we set every ring atom element to carbon
    :param input_smiles:
    :return:
    """
    mol = Chem.MolFromSmiles(input_smiles)
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)
        if atom.IsInRing():
            atom.SetAtomicNum(6)
            atom.SetIsAromatic(False)
        for bond in atom.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
    return Chem.MolToSmiles(mol)

def make_child_mol(rebuilt_smi):
    """
    Make the child molecule
    :param rebuilt_smi:
    :return:
    """
    return Chem.CanonSmiles(rebuilt_smi.replace("[Xe]","[H]"))

def get_info(atom):
    """
    Get the needed info for an atom
    :param atom:
    :return:
    """
    return [atom.GetIdx(),atom.GetNeighbors()[0].GetIdx()]


def get_type(smiles):
    """
    Get the Type of a given entity
    :param smiles: the input SMILES
    :return: the type (FG or RING)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[*;R]")):
        return "RING"
    return "FG"