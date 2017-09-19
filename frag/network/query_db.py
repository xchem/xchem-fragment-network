from neo4j.v1 import GraphDatabase
from rdkit.Chem import Draw
from rdkit import Chem
# No auth on the database
driver = GraphDatabase.driver("bolt://localhost:7687")


def find_neighbours(tx, input_str):
    return tx.run("MATCH (sta:F2 {smiles:$smiles})-[n4:F2EDGE]-(n3:F2)-[n2:F2EDGE]-(end:EM) where" \
                                                             " abs(sta.hac-end.hac) <= 3 and abs(sta.chac-end.chac) <= 1" \
                                                             " and sta.smiles <> end.smiles " \
                                                             "RETURN split(n4.label, '|')[4], split(n4.label, '|')[2], split(n2.label, '|')[2], " \
                                                             "split(n2.label, '|')[1], end.EM, end.smiles " \
                                                             "order by split(n4.label, '|')[4], split(n2.label, '|')[2];",
                  smiles=input_str)
    #return tx.run("match p = (n:F2 {smiles: $smiles})-[nm] - (m:EM) where abs(n.hac - m.hac) <= 3 and abs(n.chac - m.chac) <= 1 return split(nm.label, '|')[4], split(nm.label, '|')[1], nm.label, m.hac - n.hac, m.smiles, m.EM order by split(nm.label, '|')[4]",
    #smiles=input_str)


def get_type(r_group_form,sub_one, sub_two):
        if "." in r_group_form:
            if "C1" in sub_two:
                return "ring_linker"
            return "linker"
        if "C1" in sub_two:
            return "ring_replacement"
        return "replacement"


def write_results(out_d):
    for key in out_d:
        mols = []
        for mol in out_d[key]:
            for new_mols in out_d[key][mol]:
                mols.append(Chem.MolFromSmiles(new_mols[2]))
        # Write out the image
        img = Draw.MolsToGridImage(mols)
        img.save(key+".png")
        # out_file = open(key+".svg","w")
        # out_file.write(img)


records = []
out_d = {"replacement": {}, "ring_linker": {}, "linker": {}, "ring_replacement":{}}
with driver.session() as session:
    for record in session.read_transaction(find_neighbours,"Oc1ccc(cc1)c2ccccc2"):
        records.append(record)
        trans_from = record[1]
        trans_to = record[2]
        r_group_form = record[0]
        end_mol = record[5]
        type = get_type(r_group_form,trans_from, trans_to)
        if r_group_form not in out_d[type]:
            out_d[type][r_group_form] = [(trans_from, trans_to,end_mol)]
        else:
            out_d[type][r_group_form].append((trans_from, trans_to,end_mol))
        write_results(out_d)