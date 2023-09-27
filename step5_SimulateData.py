import os
import itertools
from functions.create_config import generate_config
import glob
import tqdm

############# USER #############
dimentions = (2048, 2048, 1024)  # in Angstroms (snowball size * snowball pixel size)
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
# generating configurations from snowballs
output_path = 'parakeet'
os.mkdir(output_path)
snowball_id = 0
for combination in itertools.product(total_doses, tilt_steps, start_angles, ice_densities, defoci):
    snowball_file = 'snowballs/{}/atomic_angposfile.txt'.format(snowball_id)
    total_dose, tilt_step, start_angle, ice_density, defocus = combination
    eletrons_per_angstrom = total_dose / len(range(start_angle, -start_angle + tilt_step, tilt_step))
    config_dir = '{}/{}'.format(output_path, snowball_id)
    os.mkdir(config_dir)
    config_path = '{}/config.yaml'.format(config_dir)
    generate_config(eletrons_per_angstrom, dimentions, phase_plates, defocus, snowball_file, start_angle, tilt_step,
                    ice_density, config_path)
    snowball_id += 1

simulation_list = glob.glob('parakeet/*')
parent_dir = os.path.abspath(os.curdir)
print(simulation_list)
# Running the simulation
for simulation in tqdm.tqdm(simulation_list):
    os.chdir(simulation)
    os.system('parakeet.run -c config.yaml --steps all')
    os.chdir(parent_dir)
