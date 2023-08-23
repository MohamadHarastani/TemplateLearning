# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Mohamad Harastani
# *
# * IGBMC, France
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'mohamad.harastani@igbmc.fr'
# *
# **************************************************************************

import mrcfile
from scipy.spatial.transform import Rotation
from scipy.ndimage import binary_dilation
import scipy
from skimage.morphology import ball
import numpy as np
from tqdm import trange
import pyfftw
import multiprocessing
import argparse
import os


def rotate3d(data, rotation, center=None, order=2):
    """Rotate a 3D volume using a rotation matrix
    @return: the data after rotation.
    """
    # Figure out the rotation center
    if center is None:
        cx = data.shape[0] / 2
        cy = data.shape[1] / 2
        cz = data.shape[2] / 2
    else:
        assert len(center) == 3
        (cx, cy, cz) = center

    from scipy import mgrid
    grid = mgrid[-cx:data.shape[0] - cx, -cy:data.shape[1] - cy, -cz:data.shape[2] - cz]
    temp = grid.reshape((3, grid.size // 3))
    temp = np.dot(np.linalg.inv(rotation), temp)
    grid = np.reshape(temp, grid.shape)
    grid[0] += cx
    grid[1] += cy
    grid[2] += cz

    # Interpolation
    from scipy.ndimage import map_coordinates
    d = map_coordinates(data, grid, order=order)
    return d


def gaussian3d(data, sigma=0):
    """Gaussian filter.

    @param data: data to be filtered.
    @param sigma: sigma of Gaussian.

    @return: filtered data.
    """
    from scipy.ndimage import gaussian_filter
    d = gaussian_filter(data, sigma)
    return d


def place(big_volume, small_volume, position):
    box_size = np.shape(small_volume)[0]
    x_start = position[0] - box_size // 2
    y_start = position[1] - box_size // 2
    z_start = position[2] - box_size // 2
    x_end = position[0] + box_size // 2
    y_end = position[1] + box_size // 2
    z_end = position[2] + box_size // 2
    big_volume[x_start:x_end, y_start:y_end, z_start:z_end] += small_volume


def unpad(big_volume, coordinates, box_size):
    # A function to keep the non-zero part of the volume, with some extra space around
    big_size = np.shape(big_volume)
    coords = np.array(coordinates)
    min_x = np.min(coords[:, 0]) - box_size
    min_y = np.min(coords[:, 1]) - box_size
    min_z = np.min(coords[:, 2]) - box_size
    max_x = np.max(coords[:, 0]) + box_size
    max_y = np.max(coords[:, 1]) + box_size
    max_z = np.max(coords[:, 2]) + box_size
    # Make sure not to cross the boundaries of the volume
    min_x = max(0, min_x)
    min_y = max(0, min_y)
    min_z = max(0, min_z)
    max_x = min(max_x, big_size[0])
    max_y = min(max_y, big_size[1])
    max_z = min(max_z, big_size[2])

    pads = [min_x, min_y, min_z, max_x, max_y, max_z]
    return pads, big_volume[min_x:max_x, min_y:max_y, min_z:max_z]


def repad(unpadded_volume, pads, actual_size):
    big_volume = np.zeros(actual_size, dtype=np.float32)
    min_x, min_y, min_z, max_x, max_y, max_z = pads
    big_volume[min_x:max_x, min_y:max_y, min_z:max_z] = unpadded_volume
    return big_volume


def snowball(mols, dim, frequencies, Iterations, coordinate_table, angle_table, output_volume, insersion_distances,
             sigma,
             gray_level_threshold, threads=None):
    """ Function to run the snowball in python
          molecule             : path for input data volume (mrc format)
          dim                  : dimensions of the big volume (phantom tomogram) in voxels [nx, ny, nz]
          Iterations           : number of molecules to attempt placing in a volume, e.g., 1000
          coordinate_table     : file path to store the output locations
          angle_table          : file path to store the output angles
          output_volume        : file path to store the output volume (phantom tomogram, mrc format)
          insersion_distances  : at what range from the bulk to insert the molecules, e.g., [5, 10]
          sigma  (double)      : sigma of the low pass filtering, e.g., 2
          gray_level_threshold : When the volumes are lowpass filtered, they are binarized. This is the threshold. e.g., 0.1
          threads              : The number of threads to use in the processing
    """
    # For debugging
    # print(mols, dim, frequencies, Iterations, coordinate_table, angle_table, output_volume, insersion_distances, sigma,
    #       gray_level_threshold, threads)

    # all coordinates are important for padding and unpadding
    all_coordinates = []

    # Monkey patch fftpack with pyfftw.interfaces.scipy_fftpack
    # Allowing the calculation of fast implementation of fft
    # Configure PyFFTW to use all cores (the default is single-threaded)
    if threads:
        pyfftw.config.NUM_THREADS = threads
    else:
        pyfftw.config.NUM_THREADS = multiprocessing.cpu_count()
    scipy.fftpack = pyfftw.interfaces.scipy_fftpack

    very_first = True

    # Get the center of the output volume and the box size of the small volume
    center = np.array(dim) // 2

    # to break out of two loops
    breaker = False
    for _ in trange(Iterations):
        if breaker:
            break
        for mol, freq in zip(mols, frequencies):
            if breaker:
                break
            # initialize angles and coordinates
            angles = []
            coordinates = []

            # read the molecule array
            molecule = []
            with mrcfile.open(mol.strip()) as mrc:
                # molecule = mrc.data
                # mrcfile reads the axis as zyx, we need xyz, so we transpose
                molecule = np.transpose(mrc.data)

            # Get the box size of the small volume
            box_size = np.shape(molecule)[0]

            if very_first:
                rotation = Rotation.random()
                # To keep the convention of relion, using ZYZ
                phi, theta, psi = rotation.inv().as_euler('ZYZ', degrees=True)
                # Rotate a copy of the molecule
                rotated_molecule = rotate3d(molecule, rotation.as_matrix())
                # Initialize a big volume
                outVolume = np.zeros(dim, dtype=np.float32)
                # place the rotated molecule in the center of the image
                place(outVolume, rotated_molecule, center)
                # Store the center and angles of the center molecule in the tables
                angles.append([phi, theta, psi])
                coordinates.append(center)
                all_coordinates.append(center)
                very_first = False

            # Actually looping:
            for _ in range(freq):
                rotation = Rotation.random()
                rotated_molecule = rotate3d(molecule, rotation.as_matrix())
                # Apply Gaussian filter to the molecule then binarize it
                rotated_binary = gaussian3d(rotated_molecule, sigma)
                rotated_binary = (rotated_binary > gray_level_threshold).astype(np.int8)
                # create the band-like shape around the molecule
                out_layer = binary_dilation(rotated_binary, ball(insersion_distances[1]))
                in_layer = binary_dilation(rotated_binary, ball(insersion_distances[0]))
                # -100 works as a penalty inside the band that should not touch any white space
                template = out_layer - 100 * in_layer
                # create a binary big volume
                # unpad the zeros and keep only the effective area
                pads, outputVol_binary = unpad(outVolume, all_coordinates, box_size)
                outputVol_binary = gaussian3d(outputVol_binary, sigma)
                outputVol_binary = (outputVol_binary > gray_level_threshold).astype(np.int8)
                # create an empty big volume that has the binary molecule in the center (template)
                # find the cross correlation between the big volume and the template
                correlation_map = scipy.signal.correlate(outputVol_binary, template, 'same', 'fft')
                correlation_map = repad(correlation_map, pads, dim)

                # set to zero the boundaries of the correlation map
                # lower bondaries
                correlation_map[:box_size // 2 + 1, :, :] = 0
                correlation_map[:, :box_size // 2 + 1, :] = 0
                correlation_map[:, :, :box_size // 2 + 1] = 0
                # higher boundaries
                correlation_map[-box_size // 2 - 1:, :, :] = 0
                correlation_map[:, -box_size // 2 - 1:, :] = 0
                correlation_map[:, :, -box_size // 2 - 1:] = 0

                # getting the coordinates of the best correlation
                coordinate = np.unravel_index(correlation_map.argmax(), correlation_map.shape)
                if np.size(coordinate) == 3:  # in case there is more than one good place to put the molecule, continue
                    if correlation_map[
                        coordinate] <= 0:  # if any part of the template overlapped with a density we stop
                        print('saturaturation is achieved')
                        breaker = True
                        break
                    else:
                        place(outVolume, rotated_molecule, coordinate)
                        coordinates.append(coordinate)
                        all_coordinates.append(coordinate)
                        angles.append(rotation.inv().as_euler('ZYZ', degrees=True))
                else:
                    continue

            # overwriting the angles and coordinates tables at every iteration
            # Saving the coordinates and angles
            angles_file = angle_table + '/{}txt'.format(os.path.basename(mol[:-3]))
            coordinates_file = coordinate_table + '/{}txt'.format(os.path.basename(mol[:-3]))
            if os.path.exists(angles_file):
                angles_prev = np.loadtxt(angles_file, skiprows=1, dtype=float, delimiter=',')
                angles = np.array(angles)
                try:
                    angles = np.concatenate((angles, angles_prev), axis=0)
                except:
                    try:
                        angles = np.vstack((angles, angles_prev))
                    except:  # keep the angles to previous
                        break

            if os.path.exists(coordinates_file):
                coordinates_prev = np.loadtxt(coordinates_file, skiprows=1, dtype=int, delimiter=',')
                try:
                    coordinates = np.concatenate((coordinates, coordinates_prev), axis=0)
                except:
                    coordinates = np.vstack((coordinates, coordinates_prev))


            np.savetxt(angles_file, angles, delimiter=',', fmt='%f', header='a_1,a_2,a_3')
            np.savetxt(coordinates_file, coordinates, delimiter=',', fmt='%d', header='c_1,c_2,c_3')

    # Saving the output volume only when the iterations are done, since we need to delete the old one
    if os.path.isfile(output_volume):
        os.system('rm {}'.format(output_volume))
    with mrcfile.new(output_volume) as mrc:
        # transposing the data to save mrc file in the right way
        # mrc.set_data(-outVolume)
        mrc.set_data(np.transpose(-outVolume))


if __name__ == '__main__':
    # this also generates --help and error handling
    CLI = argparse.ArgumentParser(
        prog='Snowball',
        description='Generates compacted distribution of input molecules in 3D'
    )
    CLI.add_argument("--mol")
    CLI.add_argument("--dim")
    CLI.add_argument("--frequencies")
    CLI.add_argument("--iterations")
    CLI.add_argument("--coordinate_table")
    CLI.add_argument("--angle_table")
    CLI.add_argument("--output_volume")
    CLI.add_argument("--insersion_distances")
    CLI.add_argument("--sigma")
    CLI.add_argument("--gray_level_threshold")
    CLI.add_argument("--threads")

    # parse the command line
    args = CLI.parse_args()
    # access CLI options
    snowball(
        args.mol[1:-1].split(','),
        np.fromstring(args.dim[1:-1], dtype=int, sep=','),
        np.fromstring(args.frequencies[1:-1], dtype=int, sep=','),
        int(args.iterations),
        args.coordinate_table,
        args.angle_table,
        args.output_volume,
        np.fromstring(args.insersion_distances[1:-1], dtype=int, sep=','),
        float(args.sigma),
        float(args.gray_level_threshold),
        None if args.threads == 'None' else int(args.threads)
    )
