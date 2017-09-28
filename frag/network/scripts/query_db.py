import argparse

from frag.utils.network_utils import get_driver,write_results
from frag.utils.vector_utils import get_exit_vector_for_xe_smi

class ReturnObject(object):

    def __init__(self, start_smi, end_smi, label, edge_count):
        """
        Build this object.
        :param start_smi:
        :param end_smi:
        :param label:
        :param frag_type:
        :param edge_count:
        """
        self.start_smi = start_smi
        self.end_smi = end_smi
        self.label = label
        self.frag_type = None
        self.edge_count = edge_count

    def __str__(self):
        out_list = [self.start_smi,self.end_smi,self.frag_type,self.edge_count,self.label]
        return str(out_list)


def find_double_edge(tx, input_str):
    return tx.run("MATCH (sta:F2 {smiles:$smiles})-[nm:F2EDGE]-(mid:F2)-[ne:F2EDGE]-(end:EM) where" \
                                                             " abs(sta.hac-end.hac) <= 3 and abs(sta.chac-end.chac) <= 1" \
                                                             " and sta.smiles <> end.smiles " \
                                                             "RETURN sta, nm, mid, ne, end " \
                                                             "order by split(nm.label, '|')[4], split(ne.label, '|')[2];",
                  smiles=input_str)

def find_proximal(tx, input_str):
    return tx.run("match p = (n:F2{smiles:$smiles})-[nm]-(m:EM) "
                  "where abs(n.hac-m.hac) <= 3 and abs(n.chac-m.chac) <= 1 "
                  "return n, nm, m "
                  "order by split(nm.label, '|')[4];",
                  smiles=input_str)


def get_type(r_group_form, sub_one, sub_two):
    if "." in r_group_form:
        if "C1" in sub_two:
            return "ring_linker"
        return "linker"
    if "C1" in sub_two:
        return "ring_replacement"
    return "replacement"


def define_double_edge_type(record):
    """
    Define the type returned for proximal systems
    :param record:
    :return:
    """
    mol_one = record["sta"]
    label = str(record["nm"]["label"].split("|")[4])
    mid_label = str(record["ne"]["label"].split("|")[4])
    mol_two = record["mid"]
    mol_three = record["end"]
    diff_one = mol_one["hac"] - mol_two["hac"]
    diff_two = mol_two["hac"] - mol_three["hac"]

    ret_obj = ReturnObject(mol_one["smiles"],mol_three["smiles"],label,2)
    if "." in label:
        ret_obj.frag_type = "LINKER"
    elif diff_one >= 0 and diff_two >=0:
        ret_obj.frag_type = "DELETION"
    elif diff_one <=0 and diff_two <=0:
        ret_obj.frag_type = "ADDITION"
    else:
        ret_obj.frag_type = "REPLACE"
    return ret_obj


def define_proximal_type(record):
    """
    Define the type returned for proximal systems
    :param record:
    :return:
    """
    mol_one = record["n"]
    label = str(record["nm"]["label"].split("|")[4])
    mol_two = record["m"]
    ret_obj = ReturnObject(mol_one["smiles"],mol_two["smiles"],label,1)
    if "." in label:
        ret_obj.frag_type = "LINKER"
    elif mol_one["hac"] - mol_two["hac"] > 0:
        ret_obj.frag_type = "DELETION"
    elif mol_one["hac"] - mol_two["hac"] < 0:
        ret_obj.frag_type = "ADDITION"
    else:
        ret_obj.frag_type = "REPLACE"
    return ret_obj

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Query a Database for a given SMILES')
    parser.add_argument('--smiles')
    args = parser.parse_args()
    driver = get_driver()
    with driver.session() as session:
        records = []
        # TODO Canonicalise inputs
        smiles = args.smiles#Chem.MolToSmiles(Chem.MolFromSmiles(args.smiles))
        for record in session.read_transaction(find_proximal,smiles):
            ans = define_proximal_type(record)
            records.append(ans)
        for record in session.read_transaction(find_double_edge, smiles):
            ans = define_double_edge_type(record)
            records.append(ans)
        for label in list(set([x.label for x in records])):
            # Linkers are meaningless
            if "." in label:
                continue
            print(label)
            print((get_exit_vector_for_xe_smi(label)))

