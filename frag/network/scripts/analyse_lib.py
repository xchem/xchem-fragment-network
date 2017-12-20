from rdkit import DataStructs
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit import Chem
from frag.network.query import get_full_graph
from frag.network.decorate import get_add_del_link


def get_sum_stats(smi_list, smiles):
    self_fp = GetMorganFingerprint(Chem.MolFromSmiles(smiles), 2)
    fps = [GetMorganFingerprint(Chem.MolFromSmiles(x),2) for x in smi_list]
    out_sims = []
    for i in fps:
        out_sims.append(DataStructs.TanimotoSimilarity(i,self_fp))
    return len(out_sims),min(out_sims),max(out_sims),sum(out_sims)/float(len(out_sims))


def run_for_smiles(smiles):
    core_list = [item for sublist in get_add_del_link(smiles) for item in sublist]
    out_dict = get_full_graph(smiles)
    new_dict = {}
    for key in out_dict:
        smi_key = key.split("_")[0]
        if smi_key in new_dict:
            new_dict[smi_key].extend(out_dict[key])
        else:
            new_dict[smi_key] = out_dict[key]

    res_dict = {}
    for smi in core_list:
        if smi in new_dict:
            res_dict[smi] = get_sum_stats(new_dict[smi],smiles)
        else:
            res_dict[smi] = None

if __name__ == "__main__":
    smiles_list = ['CN1CCN(CC1)C(=O)C2CCCN2', 'ClC=1C=CC(=CC1)C(=O)N2CCSCC2',
                   'CCC1=CC=C(S1)C(=O)NCC=2C=NN(C)C2', 'NC(=O)C1CCN(CCC=2C=CC=CC2)CC1',
                   'O=C(N1CCCCCC1)C=2C=CC=3OCOC3C2']
    res_dict = {}
    for smi in smiles_list:
        res_dict[smi] = run_for_smiles(smi)