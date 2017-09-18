from rdkit import Chem
from models import NodeHolder,Edge
from utils import make_child_mol, rebuild_smi,get_fragments


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
    # Rebuild the molecule
    rebuilt_smi = rebuild_smi(new_list,ring_ring)
    # Turn into child molecule
    child_smi = make_child_mol(rebuilt_smi)
    # Now generate the edges with input and this node
    new_node,is_new = node_holder.create_or_retrieve_node(child_smi)
    node_holder.edge_list.append(Edge(excluded_smi, rebuilt_smi, input_node, new_node))
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
    if len(fragments)<2:
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


if __name__ == "__main__":

    smiles = [x.split()[1] for x in open("../tests/data/attributes.txt").readlines()]
    smiles = ["Oc1ccc(cc1)c2ccccc2"]
    node_holder = NodeHolder()
    # Create the nodes and test with output
    for smile in smiles:
        node, is_node = node_holder.create_or_retrieve_node(smile)
        if is_node:
            create_children(node, node_holder)
    out_f = open("out_nodes.txt","w")
    for node in node_holder.node_list:
        out_f.write(str(node))
        out_f.write("\n")
    out_f = open("out_edges.txt","w")
    for edge in node_holder.get_edges():
        out_f.write(str(edge))
        out_f.write("\n")