import glob
import os
import multiprocessing
import tqdm

templates = list(glob.glob('templates/*'))
nma_paths ='templates/NMA'
os.mkdir(nma_paths)

# general parameters (keep default unless you're an expert)
modes_number = 20
cutoff = 15
spring_constant = 1
threads = multiprocessing.cpu_count()-1

for template in tqdm.tqdm(templates):
    name = os.path.basename(template)[:-4]
    nma_path = nma_paths + '/' + name
    os.mkdir(nma_path)
    args = 'prody anm {0} -s "all" --altloc "all"  --hessian --export-scipion --npzmatrices ' \
       '--npz -o {1} -p modes -n {2} -g {3} -c "{4}" -P {5} --turbo'.format(template,
                                                                    nma_path,
                                                                    modes_number,
                                                                    spring_constant,
                                                                    cutoff,
                                                                    threads)
    os.system(args)
