import os
import glob
import tqdm
import numpy
import copy

# parameters to set:
snowballs_pixel_size = 16
tomogram_pixel_size = 8
size_ratio = snowballs_pixel_size/tomogram_pixel_size
#

simulation_list = glob.glob('parakeet/*')
templates_list = [os.path.basename(template)[:-4] for template in glob.glob('templates/*.pdb')]
snowballs_path = 'snowballs/'
output_path = 'results/'

parent_dir = os.path.abspath(os.curdir)

if os.path.exists(output_path):
    print('rename or remove your previous results directory')
    exit(1)
else:
    os.mkdir(output_path)

# for every tomogram, copy the tomogram and all templates coordinates to the outpath
for simulation in tqdm.tqdm(simulation_list):
    # print(simulation)
    simulation_id = os.path.basename(simulation)
    snowball = '{}{}'.format(snowballs_path, simulation_id)
    # print(snowball)
    tomogram_path = '{}/tomogram.mrc'.format(simulation)
    # copy the tomogram path to the results:
    tomogram_output = '{}{}.mrc'.format(output_path,simulation_id)
    os.system('cp {} {}'.format(tomogram_path, tomogram_output))
    # get the coordinates ground truth
    all_coordinates_files = glob.glob('{}/coordinates/*.txt'.format(snowball))
    filenames = [os.path.basename(file)[:-4] for file in all_coordinates_files]
    templates_coordinates_files = []
    for file in all_coordinates_files:
        if os.path.basename(file)[:-4] in templates_list:
            templates_coordinates_files.append(file)
    # now loop over the template coordinates and combine them in a single file
    coordinates = None
    for file in templates_coordinates_files:
        new_coordinates = numpy.loadtxt(file, skiprows=1, delimiter=',')
        if coordinates is None:
            coordinates = copy.deepcopy(new_coordinates)
        else:
            coordinates = numpy.row_stack((coordinates,new_coordinates))
    # adjust the coordiantes to the right pixel size
    coordinates *= size_ratio
    # save the coordinates
    coordinates_output = '{}{}.txt'.format(output_path, simulation_id)
    numpy.savetxt(coordinates_output, coordinates, fmt='%f')
