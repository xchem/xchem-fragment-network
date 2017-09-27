from rdkit.Chem import AllChem
from rdkit import Chem
from frag.utils.vector_utils import get_exit_vector_for_xe_smi
from frag.utils.network_utils import get_fragments


if __name__ == "__main__":
    input_mol = Chem.MolFromSmiles("Oc1ccc(cc1)c2ccccc2")
    input_mol = AllChem.AddHs(input_mol)
    print(Chem.MolToSmiles(input_mol))
    for frag in get_fragments(input_mol):
        # TODO If it's not a ring ignore
        print frag
        print(get_exit_vector_for_xe_smi(frag))