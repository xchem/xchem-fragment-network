import argparse,random

from frag.network.query import get_picks

if __name__ == "__main__":

        parser = argparse.ArgumentParser(description='Query a Database for a given SMILES')
        parser.add_argument('--smiles')
        parser.add_argument('--num_picks')
        args = parser.parse_args()
        num_picks = int(args.num_picks)
        smiles = args.smiles
        # Add this node
        get_picks(smiles,num_picks)

