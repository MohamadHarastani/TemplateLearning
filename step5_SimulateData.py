import os
import itertools
from functions.create_config import generate_config
import tqdm
import joblib


############# USER #############
GPU_ID = ['0', '1']
dimentions = (2048, 2048, 1024)  # in Angstroms (tetris size * tetris pixel size)
phase_plates = True
# paramters used as variations (here 48 variations)
total_doses = [75, 150]
tilt_steps = [2, 4]
start_angles = [-60, -40]
ice_densities = [0.9, 1.1]
if phase_plates:
    defoci = [0, 0.5, 1]
else:
    defoci = [2.5, 3.25, 4]

############# CODE #############
# generating configurations from tetrises
output_path = 'parakeet'
os.mkdir(output_path)

# repeat the GPU_IDs so that each simulation is done on a new GPU
GPU_list = GPU_ID * 1000
tetris_id = 0

for combination in itertools.product(total_doses, tilt_steps, start_angles, ice_densities, defoci):
    tetris_file = 'tetrises/{}/atomic_angposfile.txt'.format(tetris_id)
    total_dose, tilt_step, start_angle, ice_density, defocus = combination
    eletrons_per_angstrom = total_dose / len(range(start_angle, -start_angle + tilt_step, tilt_step))
    config_dir = '{}/{}'.format(output_path, tetris_id)
    os.mkdir(config_dir)
    config_path = '{}/config.yaml'.format(config_dir)
    generate_config(eletrons_per_angstrom, dimentions, phase_plates, defocus, tetris_file, start_angle, tilt_step,
                    ice_density, config_path, GPU_list[tetris_id])
    tetris_id += 1

simulation_list = ["parakeet/{}".format(i) for i in range(tetris_id)]

parent_dir = os.path.abspath(os.curdir)
print(simulation_list)

# Running the simulation
def simulate(simulation, parent_dir):
    os.chdir(simulation)
    os.system('parakeet.run -c config.yaml --steps all')
    os.chdir(parent_dir)


joblib.Parallel(n_jobs=GPU_ID.__len__())(joblib.delayed(simulate)(simulation, parent_dir)
                                         for simulation in tqdm.tqdm(simulation_list))
