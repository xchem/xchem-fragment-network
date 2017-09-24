import argparse

from rdkit import Chem

from frag.utils.network_utils import get_driver,write_results

def find_double_edge(tx, input_str):
    return tx.run("MATCH (sta:F2 {smiles:$smiles})-[nm:F2EDGE]-(n3:F2)-[ne:F2EDGE]-(end:EM) where" \
                                                             " abs(sta.hac-end.hac) <= 3 and abs(sta.chac-end.chac) <= 1" \
                                                             " and sta.smiles <> end.smiles " \
                                                             "RETURN nm, ne, end" \
                                                             "order by split(nm.label, '|')[4], split(ne.label, '|')[2];",
                  smiles=input_str)

def find_proximal(tx, input_str):
    return tx.run("match p = (n:F2{smiles:$smiles})-[nm]-(m:EM) "
                  "where abs(n.hac-m.hac) <= 3 and abs(n.chac-m.chac) <= 1 "
                  "return nm, m "
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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Query a Database for a given SMILES')
    parser.add_argument('--smiles')
    args = parser.parse_args()
    driver = get_driver()
    with driver.session() as session:
        records = []
        out_d = {"replacement": {}, "ring_linker": {}, "linker": {}, "ring_replacement": {},"far":{}}
        for record in session.read_transaction(find_double_edge, Chem.MolToSmiles(Chem.MolFromSmiles(args.smiles))):
            trans_from = record[1]
            trans_to = record[2]
            r_group_form = record[0]
            end_mol = record[5]
            type = get_type(r_group_form,trans_from, trans_to)
            if r_group_form not in out_d[type]:
                out_d[type][r_group_form] = [(trans_from, trans_to,end_mol)]
            else:
                out_d[type][r_group_form].append((trans_from, trans_to,end_mol))
        out_d["far"] = {}
        for record in session.read_transaction(find_double_edge,Chem.MolToSmiles(Chem.MolFromSmiles(args.smiles))):
            trans_from = record[1]
            trans_to = record[2]
            r_group_form = record[0]
            end_mol = record[5]
            if r_group_form not in out_d["far"]:
                out_d["far"][r_group_form] = [(trans_from, trans_to,end_mol)]
            else:
                out_d["far"][r_group_form].append((trans_from, trans_to,end_mol))
        for record in session.read_transaction(find_proximal,Chem.MolToSmiles(Chem.MolFromSmiles(args.smiles))):
            print record
            # trans_from = record[1]
            # trans_to = record[2]
            # r_group_form = record[0]
            # end_mol = record[4]
            # if r_group_form not in out_d["near"]:
            #     out_d["near"][r_group_form] = [(trans_from, trans_to,end_mol)]
            # else:
            #     out_d["near"][r_group_form].append((trans_from, trans_to,end_mol))

        for key in out_d:
            img = write_results(out_d[key])
            for mol in img:
                out_f = open(mol+"__"+key+".svg","w")
                out_f.write(img[mol])
