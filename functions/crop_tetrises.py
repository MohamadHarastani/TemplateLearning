import glob
import os
import numpy as np
import mrcfile


def inside(new_coords, cropped_size):
    if 0 < new_coords[0] < cropped_size[0] and 0 < new_coords[1] < cropped_size[1] and 0 < new_coords[2] < cropped_size[2]:
        return True
    else:
        return False


original_size = np.array([256, 256, 64])
cropped_size = np.array([128, 128, 64])

center_shift = (original_size - cropped_size) // 2

tetrises_path = glob.glob('tetrises/*')
tetrises_path.sort()
cropped_tetrises_path = 'tetrises_cropped'

for tetris_path in tetrises_path:
    output_path = "{}/{}".format(cropped_tetrises_path, os.path.basename(tetris_path))
    # tetris_path = "tetrises/0"
    # output_path = "tetrises_cropped/0"

    if os.path.exists(output_path):
        os.system('rm -rf {}'.format(output_path))
    os.makedirs(output_path)

    input_volume = "{}/output_volume.mrc".format(tetris_path)
    output_volume = "{}/output_volume.mrc".format(output_path)

    with mrcfile.open(input_volume) as mrc:
        data = mrc.data

    data = data.transpose()
    x, y, z = center_shift // 2
    if z == 0:
        cropped_data = data[x:-x, y:-y, :]
    else:
        cropped_data = data[x:-x, y:-y, z:-z]

    cropped_data = cropped_data.transpose()
    mrcfile.write(output_volume, cropped_data)

    coordinates_path = "{}/coordinates".format(tetris_path)
    angles_path = "{}/angles".format(tetris_path)

    coordinates_files = glob.glob("{}/*.txt".format(coordinates_path))

    output_coordinates_path = "{}/coordinates".format(output_path)
    output_angles_path = "{}/angles".format(output_path)
    os.makedirs(output_coordinates_path)
    os.makedirs(output_angles_path)

    for coordinates_file in coordinates_files:
        basename = os.path.basename(coordinates_file)

        new_coordinates_file = "{}/{}".format(output_coordinates_path, basename)
        new_angles_file = "{}/{}".format(output_angles_path, basename)

        angles_file = "{}/{}".format(angles_path, basename)
        # read each the coordinates and angles file
        coords_prev = np.loadtxt(coordinates_file, skiprows=1, dtype=float, delimiter=',')
        angles_prev = np.loadtxt(angles_file, skiprows=1, dtype=float, delimiter=',')
        id = 0

        coordinates_cropped = []
        angles_cropped = []
        for coords, angles in zip(coords_prev, angles_prev):
            new_coords = coords - center_shift
            if inside(new_coords, cropped_size):
                coordinates_cropped.append(new_coords)
                angles_cropped.append(angles)
        if len(coordinates_cropped) > 0:
            np.savetxt(new_coordinates_file, coordinates_cropped, delimiter=',', fmt='%d', header='c_1,c_2,c_3')
            np.savetxt(new_angles_file, angles_cropped, delimiter=',', fmt='%f', header='a_1,a_2,a_3')
