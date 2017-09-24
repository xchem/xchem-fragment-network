import argparse
import os

from frag.utils.network_utils import get_driver

def add_nodes(tx, file_path):
    tx.run("LOAD CSV FROM '$file_path' AS line FIELDTERMINATOR "
           "' ' MERGE (:F2 { smiles: line[1], hac: toInt(line[2]), chac: toInt(line[3]), osmiles: line[4]});",
           file_path=file_path)

def add_edges(tx,file_path):
    tx.run("LOAD CSV FROM '$file_path' AS line FIELDTERMINATOR "
           "' ' MATCH (n1:F2 { smiles: line[1]}), (n2:F2 { smiles: line[2]}) MERGE (n1)-[:F2EDGE{label:line[3]}]->(n2);",
           file_path=file_path)

def add_attrs(tx,file_path):
    tx.run("LOAD CSV FROM '$file_path' AS line FIELDTERMINATOR ' ' "
           "MATCH (n:F2 { smiles: line[1]} ) set n:MOL, n:EM, n.EM=toInt(line[3]);",
           file_path=file_path)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Upload Nodes, Edges and Attributes')
    parser.add_argument('--base_dir')
    args = parser.parse_args()
    driver = get_driver()
    with driver.session() as session:
        session.write_transaction(add_nodes,os.path.join(args.base_dir,"nodes.txt"))
        session.write_transaction(add_edges,os.path.join(args.base_dir,"edges.txt"))
        session.write_transaction(add_attrs,os.path.join(args.base_dir,"attributes.txt"))
