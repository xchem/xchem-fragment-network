from neo4j.v1 import GraphDatabase

# No auth on the database
driver = GraphDatabase.driver("bolt://localhost:7687")


def find_neighbours(tx, input_str):
    return tx.run("match p = (n:F2 {smiles: $smiles})-[nm] - (m:EM) where abs(n.hac - m.hac) <= 3 and abs(n.chac - m.chac) <= 1 return split(nm.label, '|')[4], split(nm.label, '|')[1], nm.label, m.hac - n.hac, m.smiles, m.EM order by split(nm.label, '|')[4]",
    smiles=input_str)



records = []
with driver.session() as session:
    for record in session.read_transaction(find_neighbours,"Oc1ccc(cc1)c2ccccc2"):
        records.append(record)
