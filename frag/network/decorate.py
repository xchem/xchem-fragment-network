from rdkit import Chem
from rdkit.Chem import AllChem,RWMol,Atom
from frag.utils.network_utils import rebuild_smi,get_ring_ring_splits,get_fragments


def get_mol(input_smi):
    return RWMol(AllChem.AddHs(Chem.MolFromSmiles(input_smi)))

def decorate_smi(input_smi):
    """
    Decorate an input SMILES with a pseudeo molecule with all desirable changes.
    This can be added to the list of input molecules - and then we know if we get an Antimony back - it's not a real trans.
    :param input_smi: the input smiles
    :return: a list of SMILES around a given molecule.
    """
    # Add replacement groups (e.g. At) to all Ring positions in turn.
    mol = get_mol(input_smi)
    # Find matches
    patt = Chem.MolFromSmarts("[*;R]-;!@[H]")
    # Get the list of atom Indices to replace
    out_atom_repls = [x[1] for x in mol.GetSubstructMatches(patt)]
    # Now replace with At - an produce a new mol everytime
    new_mols = []
    for atom in out_atom_repls:
        rw_mol = get_mol(input_smi)
        rw_mol.ReplaceAtom(atom,Atom(85))
        newer_mol = rw_mol.GetMol()
        new_mols.append(Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(newer_mol,isomericSmiles=True)),isomericSmiles=True))
    return list(set(new_mols))


def deletion_linker_smi(input_smi):
    """
    Produce all the linker deletion SMILES and replace with Li
    :param input_smi: the input SMI
    :return:
    """
    mol = Chem.MolFromSmiles(input_smi)
    nr = mol.GetRingInfo().NumRings()
    fragments = get_fragments(mol)
    out_mols = []
    linker_mols = []
    ring_ring_splits = get_ring_ring_splits(mol)
    if ring_ring_splits:
        for ring_ring_split in ring_ring_splits:
            rebuilt_smi = rebuild_smi(ring_ring_split,ring_ring=True)
            new_mol = Chem.MolFromSmiles(rebuilt_smi)
            rebuilt_smi = rebuilt_smi.replace("Xe","Li")
            if new_mol.GetRingInfo().NumRings() < nr:
                 continue
            new_mol = link_li(rebuilt_smi)
            linker_mols.append(new_mol)
    for i in range(len(fragments)):
        new_list = []
        for j,item in enumerate(fragments):
            if i==j:
                excluded_smi = item
                continue
            new_list.append(item)
        rebuilt_smi = rebuild_smi(new_list,ring_ring=False)
        if "." in rebuilt_smi:
            new_mol = Chem.MolFromSmiles(rebuilt_smi)
            rebuilt_smi = rebuilt_smi.replace("Xe","Li")
            if new_mol.GetRingInfo().NumRings() < nr:
                continue
            new_mol = link_li(rebuilt_smi)
            linker_mols.append(new_mol)
            continue
        new_mol = Chem.MolFromSmiles(rebuilt_smi)
        # If the resulting
        if new_mol.GetRingInfo().NumRings() < nr:
            continue
        out_mols.append(new_mol)
    return (out_mols,linker_mols)

def link_li(rebuilt_smi):
    mol = Chem.MolFromSmiles(rebuilt_smi)
    mol = RWMol(mol)
    # bons =  [x[0] for x in mol.GetSubstructMatches(Chem.MolFromSmarts("[Li]"))]
    # mol.AddBond(bons[0],bons[1])
    return mol.GetMol()

def addition_smi(input_smi):
    smis = decorate_smi(input_smi)
    return [Chem.MolFromSmiles(x) for x in smis]

def get_add_del_link(smi,asSmiles=True):
    additions = addition_smi(smi)
    res = deletion_linker_smi(smi)
    linkers = res[1]
    deletions = res[0]
    if asSmiles:
        additions = [Chem.MolToSmiles(x) for x in additions]
        deletions = [Chem.MolToSmiles(x) for x in deletions]
        linkers = [Chem.MolToSmiles(x) for x in linkers]
    return [additions,deletions,linkers]