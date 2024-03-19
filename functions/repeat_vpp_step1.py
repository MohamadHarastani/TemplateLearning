# This step is not in the main pipeline, but if you have done the 10 steps with/without VPP and
# you want to use the same data to do without/with VPP
import glob
import os
import itertools
from functions.create_config import generate_config
import tqdm

# paramters to set
dimentions = (4096, 4096, 1024)  # in Angstroms
phase_plates = True

# paramters used as variations (here 48 variations)
total_doses = [75, 150]
tilt_steps = [2, 4]
start_angles = [-60, -40]
ice_densities = [0.9, 1.1]
if phase_plates:
    defoci = [0, -0.5, -1]
else:
    defoci = [-2.5, -3.25, -4]

# generating configurations from tetrises
simulation_path_path = 'parakeet'
output_path = 'vpp'
os.mkdir(output_path)
tetris_id = 0
for combination in itertools.product(total_doses, tilt_steps, start_angles, ice_densities, defoci):
    tetris_file = 'tetrises/{}/atomic_angposfile.txt'.format(tetris_id)
    total_dose, tilt_step, start_angle, ice_density, defocus = combination
    eletrons_per_angstrom = total_dose / len(range(start_angle, -start_angle + tilt_step, tilt_step))
    config_dir = '{}/{}'.format(output_path, tetris_id)
    os.mkdir(config_dir)
    config_path = '{}/config.yaml'.format(config_dir)
    generate_config(eletrons_per_angstrom, dimentions, phase_plates, defocus, tetris_file, start_angle, tilt_step,
                    ice_density, config_path)
    tetris_id += 1

# link the sample and exit wave from the previous simulations
file_ids = glob.glob('parakeet/*')
for file_id in file_ids:
    ew = '{}/exit_wave.h5'.format(file_id)
    ewl = '{}/{}/exit_wave.h5'.format(output_path,os.path.basename(file_id))
    s = '{}/sample.h5'.format(file_id)
    sl = '{}/{}/sample.h5'.format(output_path,os.path.basename(file_id))
    os.link(ew, ewl)
    os.link(s, sl)


# generate optics and image
simulation_list = glob.glob('vpp/*')
parent_dir = os.path.abspath(os.curdir)


for simulation in tqdm.tqdm(simulation_list):
    os.chdir(simulation)
    os.system('parakeet.run -c config.yaml -s sample.h5 --steps simulate.optics simulate.image')
    os.chdir(parent_dir)

# replace the old images with the new ones
# binning pixel size
PixelSize = 8  # in Angstroms, has to be integer for imod to accept it

simulation_list = glob.glob('vpp/*')
parent_dir = os.path.abspath(os.curdir)

for simulation in tqdm.tqdm(simulation_list):
    os.chdir(simulation)
    tilt_series_h5 = 'image.h5'
    tilt_series_mrc = 'image.mrc'
    if os.path.exists(tilt_series_mrc):
        os.remove(tilt_series_mrc)
    # export the data to mrc:
    os.system('parakeet.export {} -o {}'.format(tilt_series_h5, tilt_series_mrc))
    # bin the data with imod
    tilt_series_binned = 'tiltseries.mrcs'
    if os.path.exists(tilt_series_binned):
        os.remove(tilt_series_binned)
    os.system('newstack -input {} -output {} -taper 1,1 -bin {} '
              '-antialias 6 -imagebinned 1.0'.format(tilt_series_mrc, tilt_series_binned, PixelSize))
    # remove the original size tilt series
    os.system('rm {}'.format(tilt_series_mrc))
    os.chdir(parent_dir)
