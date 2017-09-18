from rdkit import Chem
from models import NodeHolder,Edge
from utils import make_child_mol, get_info


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
            mols = [x.replace("1*", "Xe") for x in Chem.MolToSmiles(nm, isomericSmiles=True).split(".")]
            out_mols.append(mols)
        return out_mols


def add_child_and_edge(new_list, input_node, excluded_smi, node_holder, ring_ring=False):
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
    new_node,is_new = node_holder.create_or_retrieve_node(child_smi)
    node_holder.edge_list.append(Edge(excluded_smi, child_smi,input_node, new_node))
    # Now generate the children for this too
    if is_new:
        create_children(new_node, node_holder)


def create_children(input_node, node_holder):
    """
    Create a series of edges from an input molecule. Iteratively
    :param input_node:
    :return:
    """
    # Get all ring-ring splits
    ring_ring_splits = get_ring_ring_splits(input_node.RDMol)
    if ring_ring_splits:
        for ring_ring_split in ring_ring_splits:
            add_child_and_edge(ring_ring_split,input_node,"[Xe]", node_holder,ring_ring=True)
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
            add_child_and_edge(new_list, input_node, excluded_smi, node_holder)



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



if __name__ == "__main__":

    node_holder = NodeHolder()
    node, is_node = node_holder.create_or_retrieve_node("Cc3sc1-n2c(C)nnc2C(CC(=O)OC(C)(C)C)N=C(c1c3C)c4ccc(Cl)cc4")
    if is_node:
        create_children(node, node_holder)
    for node in node_holder.node_list:
        print node
    for edge in node_holder.get_edges():
        print edge