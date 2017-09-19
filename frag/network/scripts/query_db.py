import json
import argparse
from frag.network.utils import get_driver,write_results
from rdkit import Chem

def find_neighbours(tx, input_str):
    return tx.run("MATCH (sta:F2 {smiles:$smiles})-[n4:F2EDGE]-(n3:F2)-[n2:F2EDGE]-(end:EM) where" \
                                                             " abs(sta.hac-end.hac) <= 3 and abs(sta.chac-end.chac) <= 1" \
                                                             " and sta.smiles <> end.smiles " \
                                                             "RETURN split(n4.label, '|')[4], split(n4.label, '|')[2], split(n2.label, '|')[2], " \
                                                             "split(n2.label, '|')[1], end.EM, end.smiles " \
                                                             "order by split(n4.label, '|')[4], split(n2.label, '|')[2];",
                  smiles=input_str)


def get_type(r_group_form, sub_one, sub_two):
    if "." in r_group_form:
        if "C1" in sub_two:
            return "ring_linker"
        return "linker"
    if "C1" in sub_two:
        return "ring_replacement"
    return "replacement"




if __name__ == "__main__":

    # Read in a SD or SMILES file - then write out into a specified directory
    parser = argparse.ArgumentParser(description='Query a Database for a given SMILES')
    parser.add_argument('--smiles')
    args = parser.parse_args()
    driver = get_driver()
    with driver.session() as session:
        records = []
        out_d = {"replacement": {}, "ring_linker": {}, "linker": {}, "ring_replacement": {}}
        for record in session.read_transaction(find_neighbours,Chem.MolToSmiles(Chem.MolFromSmiles(args.smiles))):
            trans_from = record[1]
            trans_to = record[2]
            r_group_form = record[0]
            end_mol = record[5]
            type = get_type(r_group_form,trans_from, trans_to)
            if r_group_form not in out_d[type]:
                out_d[type][r_group_form] = [(trans_from, trans_to,end_mol)]
            else:
                out_d[type][r_group_form].append((trans_from, trans_to,end_mol))
        for key in out_d:
            img = write_results(out_d[key])
            out_f = open(key+".svg","w")
            out_f.write(img)