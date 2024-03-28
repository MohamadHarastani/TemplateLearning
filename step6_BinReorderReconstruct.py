import os
import glob
import tqdm
import numpy as np
import copy
import mrcfile
from scipy.spatial.transform import Rotation
from functions.multitetris_python import place, gaussian3d, rotate3d
import argparse
from functions.utilities import createTltfile, sortTiltSeries, resize, str2bool
from termcolor import colored


def main(templates_dir, tetrises_dir, parakeet_dir, volumes_dir, output_dir, tetris_size, tetrises_pixel_size,
         tomogram_pixel_size, sigma, gray_level_threshold, adjust_to_different_size, output_size):

    tetrises_path = '{}/'.format(tetrises_dir)
    output_path = '{}/'.format(output_dir)

    output_tomograms = '{}/tomograms/'.format(output_dir)
    output_coordinates = '{}/coordinates/'.format(output_dir)
    output_segmentations = '{}/segmentations/'.format(output_dir)

    if os.path.exists(output_path):
        print(colored('rename or remove your previous results directory', 'red'))
        exit(1)
    else:
        os.mkdir(output_path)
        os.mkdir(output_tomograms)
        os.mkdir(output_coordinates)
        os.mkdir(output_segmentations)

    sirt_iterations = 0  # (advanced) zero means WBP reconstructions

    size_ratio = tetrises_pixel_size / tomogram_pixel_size
    # output tomogram thickness (be careful, this has to be equal to tetris z size at the new pixel size).
    # In other words, z = (tetris_z * tetrisPixelSize / outputPixelSize)
    z = int(size_ratio * tetris_size[2])
    PixelSize = tomogram_pixel_size

    simulation_list = glob.glob('{}/*'.format(parakeet_dir))
    parent_dir = os.path.abspath(os.curdir)

    print(colored('Binning the tilt series to the output pixel size', 'green'))
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

    parent_dir = os.path.abspath(os.curdir)

    print(colored('Sorting the tilt series before reconstruction', 'green'))
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

    print(colored('Reconstructing the tilt series into tomograms', 'green'))
    # reconstruct
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

    print(colored('Copying tomograms and coordinates to the output directory', 'green'))
    templates_list = [os.path.basename(template)[:-4] for template in glob.glob('{}/*.pdb'.format(templates_dir))]
    # for every tomogram, copy the tomogram and all templates coordinates to the outpath
    for simulation in tqdm.tqdm(simulation_list):
        # print(simulation)
        simulation_id = os.path.basename(simulation)
        tetris = '{}{}'.format(tetrises_path, simulation_id)
        # print(tetris)
        tomogram_path = '{}/tomogram.mrc'.format(simulation)
        # copy the tomogram path to the results:
        tomogram_output = '{}{}.mrc'.format(output_tomograms, simulation_id)
        os.system('cp {} {}'.format(tomogram_path, tomogram_output))
        # get the coordinates ground truth
        all_coordinates_files = glob.glob('{}/coordinates/*.txt'.format(tetris))
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
                coordinates = np.row_stack((coordinates, new_coordinates))
        # adjust the coordiantes to the right pixel size
        coordinates *= size_ratio
        # save the coordinates
        coordinates_output = '{}{}.txt'.format(output_coordinates, simulation_id)
        np.savetxt(coordinates_output, coordinates, fmt='%f')

    templates_list = [os.path.basename(template)[:-4] for template in glob.glob('{}/*.pdb'.format(templates_dir))]
    simulation_list = glob.glob('{}/*'.format(parakeet_dir))
    print(colored('Creating shape segmentation masks', 'green'))
    for simulation in tqdm.tqdm(simulation_list):
        # print(simulation)
        simulation_id = os.path.basename(simulation)
        tetris = '{}{}'.format(tetrises_path, simulation_id)
        # get the coordinates ground truth
        all_coordinates_files = glob.glob('{}/coordinates/*.txt'.format(tetris))
        all_angles_files = glob.glob('{}/angles/*.txt'.format(tetris))

        output_volume = '{}{}.mrc'.format(output_segmentations, simulation_id)
        big_volume = np.zeros(tetris_size, dtype=np.float32)
        for coordinates_file, angles_file in zip(all_coordinates_files, all_angles_files):
            basename = os.path.basename(coordinates_file)[:-4]
            if basename in templates_list:
                mol = '{}/templates/{}.mrc'.format(volumes_dir, basename)
                # read the molecule array
                molecule = []
                with mrcfile.open(mol) as mrc:
                    # mrcfile reads the axis as zyx, we need xyz, so we transpose
                    molecule = np.transpose(mrc.data)
                # reading the coordinates and angles files
                coordinates = np.atleast_2d(np.loadtxt(coordinates_file, skiprows=1, delimiter=',', dtype=int))
                angles = np.atleast_2d(np.loadtxt(angles_file, skiprows=1, delimiter=','))
                for rotation, position in zip(angles, coordinates):
                    rotation_matrix = Rotation.from_euler(angles=rotation, seq='ZYZ', degrees=True).inv().as_matrix()
                    rotated_molecule = rotate3d(molecule, rotation_matrix)
                    place(big_volume, rotated_molecule, position)

        # adjust the sampling rate
        if adjust_to_different_size:
            x, y, z = output_size
            big_volume = resize(big_volume, x, y, z)

        # binarizing the big volume
        big_volume = gaussian3d(big_volume, sigma)
        big_volume = (big_volume > gray_level_threshold).astype(np.int8)

        # save the output volume
        with mrcfile.new(output_volume) as mrc:
            # transposing the data to save mrc file in the right way
            mrc.set_data(np.transpose(big_volume))


def parse_arguments():
    parser = argparse.ArgumentParser(description='This script uses all the previous simulations to create a set of'
                                                 ' tomograms with their corresponding coordinates and segmentation'
                                                 ' masks for the templates. The output can be used to training'
                                                 ' a deep learning model (DeepFinder or similar) to pick the template'
                                                 ' in expiremental tomograms.')
    parser.add_argument('--templates', type=str, default='templates', help='Default: %(default)s. Directory where all the templates exist')
    parser.add_argument('--tetrises', type=str, default='tetrises', help='Default: %(default)s. Directory where all the tetrises exist.')
    parser.add_argument('--parakeet', type=str, default='parakeet', help='Default: %(default)s. Directory of the simulated data with parakeet.')
    parser.add_argument('--volumes', type=str, default='volumes', help='Default: %(default)s. Directory where all the volume versions of the used atomic structures exist.')
    parser.add_argument('--output', type=str, default='results', help='Default: %(default)s. Directory where the output simulations will be stored.')
    parser.add_argument('--tetris_size', nargs=3, type=int, default=[128, 128, 64],
                        help='The original tetris size in pixels (x, y, z)')
    parser.add_argument('--tetrises_pixel_size', type=int, default=16,
                        help='The original pixel size for tetris')
    parser.add_argument('--tomogram_pixel_size', type=int, default=8,
                        help='Output pixel size, in Angstroms, MUST be an integer for IMOD to accept it!')
    parser.add_argument('--output_size_is_different_than_tetris', type=str2bool, default=True,
                        help='If the pixel size of the binned tomograms is different from the tetris (which is usually the case)')
    parser.add_argument('--output_size', nargs=3, type=int, default=[256, 256, 128],
                        help='Output tomogram size after binning (size of tetris in Angstroms over the binning factor) in pixels (x, y, z)')
    parser.add_argument('--sigma', type=float, default=0.2,
                        help='Sigma value for Gaussian filter (used during binarization, the bigger it is the bigger the shape mask)')
    parser.add_argument('--gray_level_threshold', type=int, default=100,
                        help='Gray level threshold for segmentation (the higher it is, the more blob-ish the segementation.)'
                             ' Usually, keeping it close to zero gives a good segmentation. If the model does not train, try incresing it.')

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    main(args.templates, args.tetrises, args.parakeet, args.volumes, args.output, args.tetris_size, args.tetrises_pixel_size,
         args.tomogram_pixel_size, args.sigma, args.gray_level_threshold, args.output_size_is_different_than_tetris,
         args.output_size)