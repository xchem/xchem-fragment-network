import argparse
import os

from rdkit import Chem

from frag.network.models import NodeHolder,Edge,Attr
from frag.utils.network_utils import build_network, write_data
from tqdm import tqdm
if __name__ == "__main__":

    # Read in a SD or SMILES file - then write out into a specified directory
    parser = argparse.ArgumentParser(description='Convert a SMILES or SDFile to input for Astex Fragment network.')
    parser.add_argument('--input')
    parser.add_argument('--num_splits')
    args = parser.parse_args()
    num_splits = int(args.num_splits)
    out_files = []
    for i in range(num_splits):
        out_files.appned(Chem.SDWriter(args.input.replace(".sdf","_"+str(i+1)+".sdf")))
    attrs = []
    counter =0
    for x in tqdm(Chem.SDMolSupplier(args.input)):
        if x is None:
            continue
        id_this = counter % num_splits
        out_files[id_this].write(x)
        counter+=1
