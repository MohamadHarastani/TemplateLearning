import os
import glob
import tqdm
import yaml
import numpy as np
import copy
import mrcfile
from scipy.spatial.transform import Rotation
from functions.multisnowball_python import place, gaussian3d, rotate3d

############# USER #############
snowball_size = (128, 128, 64)
snowballs_pixel_size = 16
tomogram_pixel_size = 8  # output pixel size, in Angstroms, has to be integer for imod to accept it
sirt_iterations = 0  # set to zero if you want only WBP reconstructions
# durning the generations of segmentations (use the same parameters as snowball)
sigma = 2
gray_level_threshold = 100
# Match the output segmentation to the tomogram pixel size
adjust_to_different_size = True
output_size = (256, 256, 128)

############# CODE #############
size_ratio = snowballs_pixel_size/tomogram_pixel_size
z= int(size_ratio*snowball_size[2])  # output tomogram thickness (be careful, this has to be equal to snowball z size at the new pixel size). In other words, z = (snowball_z * snowballPixelSize / outputPixelSize)
PixelSize = tomogram_pixel_size

simulation_list = glob.glob('parakeet/*')
parent_dir = os.path.abspath(os.curdir)

print(simulation_list)

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


def createTltfile(name):
    presorted = '{}/presorted.tlt'.format(name)
    # read the config and get starting and step
    config_file = '{}/config.yaml'.format(name)
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    start_angle = int(config['scan']['start_angle'])
    step_angle = int(config['scan']['step_angle'])
    # create and angles list
    angles = []
    for i in range(0, -start_angle + step_angle, step_angle):
        if i == 0:
            angles.append(i)
        else:
            angles.append(-i)
            angles.append(i)
    # save as a tlt file
    if os.path.exists(presorted):
        os.remove(presorted)
    np.savetxt(presorted, angles, fmt='%f')
    return presorted


def sortTiltSeries(tsfn, tltfn, output):
    # Unstacking the tilt series
    args = 'e2spt_tiltstacker.py --unstack {} --tltfile {}'.format(tsfn, tltfn)
    os.system(args)
    tilts = np.loadtxt(tltfn)
    framesfn = sorted(list(glob.glob('sptstacker_01/*.mrcs')))
    indices = sorted(range(len(tilts)), key=lambda index: tilts[index])
    sorted_tltfile = [list(tilts)[index] for index in indices]
    sorted_fns = [framesfn[index] for index in indices]

    # now creating an imod file:
    np.savetxt('sorted.tlt', np.array(sorted_tltfile), fmt="%.2f")
    with open('imod_file.txt', 'x') as f:
        f.write('{}\n'.format(len(tilts)))
        for fn in sorted_fns:
            f.write(fn + '\n')
            f.write('0\n')

    # removing the old link to the unsorted tilt series and creating the sorted tilt series stack
    command = 'newstack -fileinlist {} -output {}'.format('imod_file.txt', output)
    os.system(command)


simulation_list = glob.glob('parakeet/*')
parent_dir = os.path.abspath(os.curdir)

print(simulation_list)

# sort tilt series
for simulation in tqdm.tqdm(simulation_list):
    presorted = createTltfile(simulation)
    os.chdir(simulation)
    tsfn = 'tiltseries.mrcs'
    output = 'sorted.mrcs'
    # if the file exists, delete it
    if os.path.exists('imod_file.txt'):
        os.system('rm -rf {} {}'.format('imod_file.txt', 'sptstacker*'))
    sortTiltSeries(tsfn, os.path.basename(presorted), output)
    os.chdir(parent_dir)

# # reconstruct
for simulation in tqdm.tqdm(simulation_list):
    os.chdir(simulation)
    tsfn = 'sorted.mrcs'
    reconstruction = 'sorted.rec'
    tomogram = 'tomogram.mrc'
    tltfile = 'sorted.tlt'
    args = 'tilt -InputProjections {} -OutputFile {} -TILTFILE {} -THICKNESS {}' \
           ' -FalloffIsTrueSigma 1 -RADIAL 0.5,0.5 -SHIFT 0.0,0.0 -OFFSET 0.0,0.0 -MODE 1 -PERPENDICULAR' \
           ' -AdjustOrigin -UseGPU 0 -ActionIfGPUFails 2,2'.format(tsfn, reconstruction, tltfile, z)
    if sirt_iterations != 0:
        args += ' -FakeSIRTiterations {}'.format(sirt_iterations)
    os.system(args)
    args = 'trimvol -rx {} {}'.format(reconstruction, tomogram)
    os.system(args)
    os.remove(reconstruction)
    os.chdir(parent_dir)


simulation_list = glob.glob('parakeet/*')
templates_list = [os.path.basename(template)[:-4] for template in glob.glob('templates/*.pdb')]
snowballs_path = 'snowballs/'
output_path = 'results/'
output_tomograms = 'results/tomograms/'
output_coordinates = 'results/coordinates/'
output_segmentations = 'results/segmentations/'

parent_dir = os.path.abspath(os.curdir)

if os.path.exists(output_path):
    print('rename or remove your previous results directory')
    exit(1)
else:
    os.mkdir(output_path)
    os.mkdir(output_tomograms)
    os.mkdir(output_coordinates)
    os.mkdir(output_segmentations)

# for every tomogram, copy the tomogram and all templates coordinates to the outpath
for simulation in tqdm.tqdm(simulation_list):
    # print(simulation)
    simulation_id = os.path.basename(simulation)
    snowball = '{}{}'.format(snowballs_path, simulation_id)
    # print(snowball)
    tomogram_path = '{}/tomogram.mrc'.format(simulation)
    # copy the tomogram path to the results:
    tomogram_output = '{}{}.mrc'.format(output_tomograms,simulation_id)
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
        new_coordinates = np.loadtxt(file, skiprows=1, delimiter=',')
        if coordinates is None:
            coordinates = copy.deepcopy(new_coordinates)
        else:
            coordinates = np.row_stack((coordinates,new_coordinates))
    # adjust the coordiantes to the right pixel size
    coordinates *= size_ratio
    # save the coordinates
    coordinates_output = '{}{}.txt'.format(output_coordinates, simulation_id)
    np.savetxt(coordinates_output, coordinates, fmt='%f')


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

    output_volume = '{}{}.mrc'.format(output_segmentations, simulation_id)
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
