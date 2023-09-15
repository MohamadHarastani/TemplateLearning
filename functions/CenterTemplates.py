import glob
import numpy as np
import tqdm
import prody

templates = list(glob.glob('templates/*.cif'))

for template in tqdm.tqdm(templates):
    pdb = prody.parsePDB(template)
    pdb._coords -= np.mean(pdb._coords, axis=1)
    new_name = template[:-4] + '_centered'
    prody.writePDB(new_name, pdb)
