from rdkit import Chem
from neo4j.v1 import GraphDatabase
from rdkit.Chem import MCS,AllChem,Draw


SMARTS_PATTERN = "[*;R]-;!@[*]"

def write_results(input_dict):
    """
    Write out a structure list of results.
    :param input_dict:
    :return:
    """
    mols = []
    out_imgs = {}
    for mol in input_dict:
        out_imgs[mol] = []
        for new_mols in input_dict[mol]:
            m = Chem.MolFromSmiles(new_mols[2])
            mols.append(m)
        if len(mols) > 2:
            p = Chem.MolFromSmarts(MCS.FindMCS(mols).smarts)
            AllChem.Compute2DCoords(p)
            for m in mols: AllChem.GenerateDepictionMatching2DStructure(m, p)
        out_imgs[mol] = Draw.MolsToGridImage(mols,useSVG=True)
        # Write out the image
    return out_imgs


def get_fragments(input_mol):
    """
    Find the frgments for a given molecule
    :param input_mol:
    :return:
    """
    atom_indices = input_mol.GetSubstructMatches(Chem.MolFromSmarts(SMARTS_PATTERN))
    if  atom_indices:
        counter = 100
        labels = []
        bs = []
        for bi in atom_indices:
            b = input_mol.GetBondBetweenAtoms(bi[0],bi[1])
            labels.append((counter,counter))
            bs.append(b.GetIdx())
            counter += 1
        input_mol = Chem.FragmentOnBonds(input_mol,bs,dummyLabels=labels)
    return [x.replace("*","Xe") for x in Chem.MolToSmiles(input_mol,isomericSmiles=True).split(".")]


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
    return Chem.MolToSmiles(mol,isomericSmiles=True)


def recombine_edges(output_edges):
    """
    Recombine a list of edges based on their rules.
    Recombines identical Xe isotopes. Remove isotopes.
    :param output_edges:
    :return:
    """
    mol = Chem.MolFromSmiles(".".join(output_edges))
    # Dictionary of atom's to bond together and delete if they come in pairs
    iso_dict = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==54:
            # Get the isotope
            iso = atom.GetIsotope()
            if iso in iso_dict:
                iso_dict[iso].append(get_info(atom))
            else:
                iso_dict[iso] = [get_info(atom)]
    mw = Chem.RWMol(mol)
    # Add bonds first
    del_indices = []
    for isotope in iso_dict:
        if len(iso_dict[isotope])>1:
            mw.AddBond(iso_dict[isotope][0][1],
                       iso_dict[isotope][1][1],
                       Chem.BondType.SINGLE)
            del_indices.append(iso_dict[isotope][0][0])
            del_indices.append(iso_dict[isotope][1][0])
    # Now delete atoms
    del_count = 0
    for atom_index in sorted(del_indices):
        mw.RemoveAtom(atom_index-del_count)
        del_count+=1
    Chem.SanitizeMol(mw)
    return Chem.MolToSmiles(mw,isomericSmiles=True)


def rebuild_smi(input_list, ring_ring):
    """
    Rebuild a SMILES
    :param input_list: the list of fragments to be rebuilt
    :param ring_ring: a boolean - indicating if this is a ring ring split
    :return:
    """
    if ring_ring:
        rebuilt_smi = ".".join(input_list)
    else:
        rebuilt_smi = recombine_edges(input_list)
    return rebuilt_smi


def make_child_mol(rebuilt_smi):
    """
    Make the child molecule
    :param rebuilt_smi:
    :return:
    """
    return Chem.CanonSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(rebuilt_smi)).replace("[Xe]","[H]"))


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


def get_driver():
    """
    Get the driver to the network connection
    :return: the driver for the graphdabase
    """
    # No auth on the database
    driver = GraphDatabase.driver("bolt://localhost:7687")
    return driver
