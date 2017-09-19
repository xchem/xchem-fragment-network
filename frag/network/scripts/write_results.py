from rdkit.Chem import Draw
from rdkit import Chem
import json

def write_results(out_d):
    for key in out_d:
        mols = []
        for mol in out_d[key]:
            for new_mols in out_d[key][mol]:
                mols.append(Chem.MolFromSmiles(new_mols[2]))
        # Write out the image
        img = Draw.MolsToGridImage(mols,useSVG=True)
        # img.save(key + ".png")
        out_file = open(key+".svg","w")
        out_file.write(img)

if __name__ == "__main__":
    write_results(json.load(open("data.json")))

