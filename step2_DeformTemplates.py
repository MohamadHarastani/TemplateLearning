import glob
import os
import random
import tqdm
import numpy
import prody

############# USER #############
# this range is good for a template of size 200 kDa, and it can be good for any template. If you want more variations,
# increase this number (e.g. double it)
deformation_range = 100
deformations_per_template = 3

############# CODE #############
templates = list(glob.glob('templates/*.pdb'))
nma_paths ='templates/NMA/'

for template in tqdm.tqdm(templates):
    name = os.path.basename(template)[:-4]
    nma_vectors_path = glob.glob(nma_paths + name + '/modes/*')
    modes = []
    for vector in nma_vectors_path:
        # ignore the first six modes that correspond to rigid-body deformations and read the other modes
        if vector.endswith(('.1', '.2', '.3', '.4', '.5', '.6')):
            continue
        else:
            mode = numpy.loadtxt(vector)
        modes.append(mode)
    # now use random combinations of the modes to generate a 1000 deformations

    for i in range(deformations_per_template):
        pdb = prody.parsePDB(template)
        pdb_ori = prody.parsePDB(template)
        for mode in modes:
            mode = numpy.array(mode)
            # deform using the mode in a random range
            amplitude = (random.random() - 0.5) * 2 * deformation_range
            mode *= amplitude
            pdb._coords += mode
        new_name = 'templates/{}_{}.pdb'.format(name, i)
        prody.writePDB(new_name, pdb)
