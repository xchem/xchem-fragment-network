from rdkit import Chem
from rdkit.Chem import AllChem,RWMol,Atom
from frag.utils.network_utils import SMARTS_PATTERN
from rdkit.Chem.Scaffolds import MurckoScaffold


def get_mol(input_smi):
    return RWMol(AllChem.AddHs(Chem.MolFromSmiles(input_smi)))

def decorate_smi(input_smi):
    """
    Decorate an input SMILES with a pseudeo molecule with all desirable changes.
    This can be added to the list of input molecules - but we know if we get an Antimony back - it's not a real trans.
    :param input_smi:
    :return:
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

