import os
import glob
import tqdm

simulation_list = glob.glob('parakeet/*')
parent_dir = os.path.abspath(os.curdir)

print(simulation_list)

for simulation in tqdm.tqdm(simulation_list):
    os.chdir(simulation)
    os.system('parakeet.run -c config.yaml --steps all')
    os.chdir(parent_dir)
