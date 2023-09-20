import os
import glob
import tqdm
import numpy as np
import copy
import mrcfile
from scipy.spatial.transform import Rotation

from functions.multisnowball_python import place, gaussian3d, rotate3d

# parameters to set:
snowball_size = (256, 256, 64)
sigma = 2
gray_level_threshold = 100
#
adjust_to_different_size = True
output_size = (316, 316, 79)


def resize(data, x, y, z):
    """Resize the data.

    @param data: input data.
    @param x: resized dimension x.
    @param y: resized dimension y.
    @param z: resized dimension z.

    @return: resized data.
    """
    s = data.shape
    from scipy import mgrid, array
    from scipy.ndimage import map_coordinates
    grid = mgrid[0:s[0] - 1:x * 1j, 0:s[1] - 1:y * 1j, 0:s[2] - 1:z * 1j]
    d = map_coordinates(data, grid, order=0)
    return d


simulation_list = glob.glob('parakeet/*')
templates_list = [os.path.basename(template)[:-4] for template in glob.glob('templates/*.pdb')]
snowballs_path = 'snowballs/'
output_path = 'segmentations/'

parent_dir = os.path.abspath(os.curdir)

if os.path.exists(output_path):
    print('rename or remove your previous segmentation directory')
    exit(1)
else:
    os.mkdir(output_path)

# for every tomogram, copy the tomogram and all templates coordinates to the outpath
for simulation in tqdm.tqdm(simulation_list):
    # print(simulation)
    simulation_id = os.path.basename(simulation)
    snowball = '{}{}'.format(snowballs_path, simulation_id)
    # print(snowball)
    # copy the tomogram path to the results:
    # get the coordinates ground truth
    all_coordinates_files = glob.glob('{}/coordinates/*.txt'.format(snowball))
    all_angles_files = glob.glob('{}/angles/*.txt'.format(snowball))

    coordinates_filenames = [os.path.basename(file)[:-4] for file in all_coordinates_files]
    angles_filenames = [os.path.basename(file)[:-4] for file in all_angles_files]

    output_volume = '{}/{}.mrc'.format(output_path, simulation_id)
    big_volume = np.zeros(snowball_size, dtype=np.float32)
    for coordinates_file, angles_file in zip(all_coordinates_files, all_angles_files):
        basename = os.path.basename(coordinates_file)[:-4]
        if basename in templates_list:
            # print(coordinates_file)
            # print(angles_file)
            # print('volumes/templates/{}.mrc'.format(basename))
            mol = 'volumes/templates/{}.mrc'.format(basename)
            # read the molecule array
            molecule = []
            with mrcfile.open(mol) as mrc:
                # molecule = mrc.data
                # mrcfile reads the axis as zyx, we need xyz, so we transpose
                molecule = np.transpose(mrc.data)

            # reading the coordinates and angles files
            coordinates = np.loadtxt(coordinates_file, skiprows=1, delimiter=',', dtype=int)
            angles = np.loadtxt(angles_file, skiprows=1, delimiter=',')
            for rotation, position in zip(angles, coordinates):
                rotation_matrix = Rotation.from_euler(angles=rotation, seq='ZYZ', degrees=True).inv().as_matrix()
                rotated_molecule = rotate3d(molecule, rotation_matrix)
                # Apply Gaussian filter to the molecule then binarize it
                rotated_binary = gaussian3d(rotated_molecule, sigma)
                rotated_binary = (rotated_binary > gray_level_threshold).astype(np.int8)
                place(big_volume, rotated_molecule, position)

    # adjust the sampling rate
    if adjust_to_different_size:
        x, y, z = output_size
        big_volume = resize(big_volume, x, y, z)
        big_volume = (big_volume > gray_level_threshold).astype(np.int8)
    # save the output volume

    with mrcfile.new(output_volume) as mrc:
        # transposing the data to save mrc file in the right way
        # mrc.set_data(-outVolume)
        mrc.set_data(np.transpose(big_volume))
