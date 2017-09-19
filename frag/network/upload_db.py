from neo4j.v1 import GraphDatabase

# No auth on the database
driver = GraphDatabase.driver("bolt://localhost:7687")

def add_node(tx, smiles,hac,chac,osmiles):
    tx.run("MERGE (:F2 { smiles: $smiles, hac: toInt($hac), chac: toInt($chac), osmiles: $osmiles})",
           smiles=smiles,hac=hac,chac=chac,osmiles=osmiles)

def add_edge(tx,smiles,smiles_two,edge_meta):
    tx.run("MATCH (n1:F2 { smiles: $smiles}), (n2:F2 { smiles: $smiles_two}) MERGE (n1)-[:F2EDGE{label:$edge_meta}]->(n2)",
           smiles=smiles, smiles_two=smiles_two, edge_meta=edge_meta)

def add_attr(tx,smiles,attr):
    tx.run("MATCH (n:F2 { smiles: $smiles} ) set n:MOL, n:EM, n.EM=$attr",
           smiles=smiles,attr=attr)


with driver.session() as session:
    for line in open("../tests/data/nodes.txt").readlines():
        session.write_transaction(add_node,line.split()[1],line.split()[2],line.split()[3],line.split()[4])
    for line in open("../tests/data/edges.txt").readlines():
        session.write_transaction(add_edge,line.split()[1],line.split()[2],line.split()[3])
    for line in open("../tests/data/attributes.txt").readlines():
        session.write_transaction(add_attr,line.split()[1],line.split()[3])