from rdkit import Chem

node_list = []
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


def get_fragments(input_mol):
    """
    Find the frgments for a given molecule
    :param input_mol:
    :return:
    """
    smarts_pattern = "[*;R]-;!@[*]"
    bond_indices = input_mol.GetSubstructMatches(Chem.MolFromSmarts(smarts_pattern))
    if not bond_indices:
        return None
    counter = 100
    labels = []
    bs = []
    for bi in bond_indices:
        b = input_mol.GetBondBetweenAtoms(bi[0],bi[1])
        labels.append((counter,counter))
        bs.append(b.GetIdx())
        counter += 1
    nm = Chem.FragmentOnBonds(input_mol,bs,dummyLabels=labels)
    return [x.replace("*","Xe") for x in Chem.MolToSmiles(nm,isomericSmiles=True).split(".")]


def get_ring_ring_splits(input_mol):
    """
    Get and break Atom-Atom pairs in two different rings.
    :param input_mol:
    :return:
    """
    #TODO Fix for fused e.g.s
    RI = input_mol.GetRingInfo()
    rings = RI.AtomRings()
    out_mols = []
    bonds = [item for sublist in RI.BondRings() for item in sublist]
    bs = []
    for bond in input_mol.GetBonds():
        if bond.GetIdx() in bonds: continue
        id_one = bond.GetBeginAtomIdx()
        id_two = bond.GetEndAtomIdx()
        # Now find all pairs that are in both
        for ring in rings:
            if id_one in ring:
                for ring_two in rings:
                    if ring==ring_two:
                        continue
                    if id_two in ring_two:
                        bs.append(bond.GetIdx())
    if bs:
        for b in bs:
            nm = Chem.FragmentOnBonds(input_mol,[b],dummyLabels=[(1,1)])
            # Only takes first
            # TODO Take all possibilities
            mols = [x.replace("1*", "Xe") for x in Chem.MolToSmiles(nm, isomericSmiles=True).split(".")]
            out_mols.append(mols)
        return out_mols


def make_child_mol(rebuilt_smi):
    """
    Make the child molecule
    :param rebuilt_smi:
    :return:
    """
    return Chem.CanonSmiles(rebuilt_smi.replace("[Xe]","[H]"))


def add_child_and_edge(new_list, input_node, excluded_smi,ring_ring=False):
    """
    :param input_pair:
    :return:
    """
    if ring_ring:
        rebuilt_smi = ".".join(new_list)
    else:
        rebuilt_smi = recombine_edges(new_list)
    child_smi = make_child_mol(rebuilt_smi)
    # Now generate the edges with input and this node
    new_edge = Edge(excluded_smi, make_child_mol(rebuilt_smi))
    input_node.EDGES.append(new_edge)
    new_node,is_new = create_or_retrieve_node(child_smi)
    new_node.EDGES.append(new_edge)
    new_edge.NODES = [input_node, new_node]
    # Now generate the children for this too
    if is_new:
        create_children(new_node)


def create_or_retrieve_node(child_smi):
    new_node = Node(Chem.MolFromSmiles(child_smi))
    new_list = [x for x in node_list if x == new_node]
    if new_list:
        return new_list[0],False
    node_list.append(new_node)
    return new_node,True


def create_children(input_node):
    """
    Create a series of edges from an input molecule. Iteratively
    :param input_node:
    :return:
    """
    # Get all ring-ring splits
    ring_ring_splits = get_ring_ring_splits(input_node.RDMol)
    if ring_ring_splits:
        for ring_ring_split in ring_ring_splits:
            add_child_and_edge(ring_ring_split,input_node,"[Xe]",ring_ring=True)
    fragments = get_fragments(input_node.RDMol)
    if not fragments:
        return
    # Now remove one item on each iteration
    for i in range(len(fragments)):
        new_list = []
        for j,item in enumerate(fragments):
            if i==j:
                excluded_smi = item
                continue
            new_list.append(item)
            add_child_and_edge(new_list, input_node, excluded_smi)



def get_info(atom):
    """
    Get the needed info for an atom
    :param atom:
    :return:
    """
    return [atom.GetIdx(),atom.GetNeighbors()[0].GetIdx()]


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
    return Chem.MolToSmiles(mw)


class Node(object):
    """
    A Class to hold a Node on the graph
    """
    def __eq__(self, other):
        return self.SMILES==other.SMILES

    def __init__(self,input_mol):
        if type(input_mol) == str:
            input_mol = Chem.MolFromSmiles(input_mol)
        self.SMILES = Chem.MolToSmiles(input_mol, isomericSmiles=True)
        self.HAC = input_mol.GetNumHeavyAtoms()
        self.RAC,split_indices = get_num_ring_atoms(input_mol)
        self.RING_SMILES = simplified_graph(self.SMILES)
        self.RDMol = input_mol
        self.EDGES = []


def get_type(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[*;R]")):
        return "RING"
    return "FG"

class Edge(object):
    """
    A Class to hold an edge
    """
    def __init__(self, excluded_smi,rebuilt_smi):
        self.EXCLUDE_SMILES = excluded_smi
        self.EXCLUDE_TYPE = get_type(excluded_smi)
        self.REBUILT_SMILES = rebuilt_smi
        self.REBUILT_RING_SMILES = simplified_graph(rebuilt_smi)
        self.REBUILT_TYPE = get_type(rebuilt_smi)
        self.EXCLUDED_RING_SMILES = simplified_graph(excluded_smi)
        self.NODES = []




if __name__ == "__main__":

    node = Node("Oc1ccc(cc1)c2ccccc2c3ccccc3")
    node_list.append(node)
    create_children(node)
    for node in node_list:
        print "NODE", node.SMILES, str(node.HAC), str(node.RAC), node.RING_SMILES
        for edge in node.EDGES:
            print "EDGE",edge.NODES[0].SMILES,edge.NODES[1].SMILES,\
                edge.EXCLUDE_TYPE, edge.EXCLUDE_SMILES,edge.EXCLUDED_RING_SMILES,\
                edge.REBUILT_TYPE,edge.REBUILT_SMILES,edge.REBUILT_RING_SMILES