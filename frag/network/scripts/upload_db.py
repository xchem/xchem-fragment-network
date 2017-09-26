import argparse
import os

from frag.utils.network_utils import get_driver
from tqdm import tqdm

def add_node(tx, smiles,hac,chac,osmiles):
    tx.run("MERGE (:F2 { smiles: $smiles, hac: toInt($hac), chac: toInt($chac), osmiles: $osmiles})",
           smiles=smiles,hac=hac,chac=chac,osmiles=osmiles)

def add_edge(tx,smiles,smiles_two,edge_meta):
    tx.run("MATCH (n1:F2 { smiles: $smiles}), (n2:F2 { smiles: $smiles_two}) MERGE (n1)-[:F2EDGE{label:$edge_meta}]->(n2)",
           smiles=smiles, smiles_two=smiles_two, edge_meta=edge_meta)

def add_attr(tx,smiles,attr):
    tx.run("MATCH (n:F2 { smiles: $smiles} ) set n:MOL, n:EM, n.EM=$attr",
           smiles=smiles,attr=attr)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Upload Nodes, Edges and Attributes')
    parser.add_argument('--base_dir')
    args = parser.parse_args()
    driver = get_driver()

    with driver.session() as session:

        for line in tqdm(open(os.path.join(args.base_dir,"nodes.txt")).readlines()):
            session.write_transaction(add_node,line.split()[1],line.split()[2],line.split()[3],line.split()[4])
        for line in tqdm(open(os.path.join(args.base_dir,"edges.txt")).readlines()):
            session.write_transaction(add_edge,line.split()[1],line.split()[2],line.split()[3])
        for line in tqdm(open(os.path.join(args.base_dir,"attributes.txt")).readlines()):
            session.write_transaction(add_attr,line.split()[1],line.split()[3])