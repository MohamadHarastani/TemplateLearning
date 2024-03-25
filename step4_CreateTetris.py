import os
import glob
import tqdm
import numpy
import random
from scipy.spatial.transform import Rotation
import argparse
from functions.utilities import str2bool

def main(volumes_dir, templates_dir, distractors_dir, output, frequencies_csv, tetris_sampling_rate, number_of_tetrises, dimensions,
         density_ratio, iterations, insertion_distances, sigma, gray_level_threshold, threads, grind):

    templates = list(glob.glob('{}/templates/*.mrc'.format(volumes_dir)))
    # repleat the templates just in case a template is alone, or the number of templates are less than distractors
    templates = templates * 200

    distractors_and_frequencies = numpy.loadtxt(frequencies_csv, delimiter=',', dtype='str')
    distractors_and_frequencies = distractors_and_frequencies.transpose()

    if not os.path.exists(output):
        os.mkdir(output)

    for tetris_number in tqdm.tqdm(range(number_of_tetrises)):
        molecules_list = []
        tetris_path = '{}/{}'.format(output, tetris_number)
        coordinates_tables_path = '{}/coordinates'.format(tetris_path)
        angles_tables_path = '{}/angles'.format(tetris_path)
        output_volume_path = '{}/output_volume.mrc'.format(tetris_path)
        os.mkdir(tetris_path)
        os.mkdir(coordinates_tables_path)
        os.mkdir(angles_tables_path)

        random.shuffle(templates)
        # fix the names to include the path and _centered.mrc, and add the templates
        for i in range(len(distractors_and_frequencies)):
            new_line = ['{}/distractors/{}_centered.mrc'.format(volumes_dir, distractors_and_frequencies[i, 0]),
                        distractors_and_frequencies[i, 1]]
            molecules_list.append(new_line)
            if i % density_ratio == 0:
                new_line = [templates[i], distractors_and_frequencies[i, 1]]

            molecules_list.append(new_line)

        # now run the tetris
        molecules = list(numpy.array(molecules_list)[:, 0])
        frequencies = list(numpy.array(molecules_list)[:, 1])
        # TODO: make sure the function is in the path
        arg = "python functions/multitetris_python.py --mol '{}' --dim '{}' --frequencies '{}' --iterations {}" \
              " --coordinate_table {} --angle_table {} --output_volume {} --insersion_distances '{}' --sigma {}" \
              " --gray_level_threshold {} --threads {} --grind {}".format(molecules, dimensions, frequencies, iterations,
                                                               coordinates_tables_path, angles_tables_path,
                                                               output_volume_path, insertion_distances, sigma,
                                                               gray_level_threshold, threads, grind)
        os.system(arg)

    print('Done with generating tetrises, now adapting the file convention to Parakeet')
    # converting the tetrises coordinates to parakeet coordinates
    tetrises_path = glob.glob('{}/*'.format(output))
    templates = glob.glob('{}/*'.format(templates_dir))
    distractors = glob.glob('{}/*'.format(distractors_dir))

    templates_base_names = [os.path.basename(template)[:-4] for template in templates]
    distractors_base_names = [os.path.basename(distractor)[:-4] for distractor in distractors]


    for tetris_path in tqdm.tqdm(tetrises_path):
        coordinates_tables = glob.glob('{}/coordinates/*.txt'.format(tetris_path))
        angles_tables = glob.glob('{}/angles/*.txt'.format(tetris_path))
        txt = ''
        for coordinates_table, angles_tables in zip(coordinates_tables, angles_tables):
            basename = os.path.basename(coordinates_table)[:-4]
            if basename in templates_base_names:
                filename = '{}/{}.pdb'.format(os.path.abspath(templates_dir), basename)
            elif basename in distractors_base_names:
                filename = '{}/{}.pdb'.format(os.path.abspath(distractors_dir), basename)
            else:
                exit()
            """ Warning! difficult text file handling """
            # [{position: [1, 2, 3], orientation: [4, 5, 6]}, {position: [4, 4, 5], orientation: [4, 7, 7]}, ...]

            # Reading the lines of coordinates and angles, ignoring the first and last lines (first has the letters and last is
            # empty)
            coorinates = open(coordinates_table).read().split('\n')[1:-1]
            angles = open(angles_tables).read().split('\n')[1:-1]

            txt += '      - filename: {}\n'.format(filename)
            txt += '        instances: ['
            for coords, angs in zip(coorinates, angles):
                # multiplying the coordinates by the voxel size
                coords = [float(i) * tetris_sampling_rate for i in coords.split(',')]
                txt += '{position: [' + '{}'.format(coords)[1:-1]
                # converting the angles to rotation vectors
                angs = [float(i) for i in angs.split(',')]
                T = Rotation.from_euler('ZYZ', [angs[0], angs[1], angs[2]], degrees=True).inv().as_matrix()
                angs = list(Rotation.from_matrix(T).as_rotvec())
                txt += '], orientation: [' + '{}'.format(angs)[1:-1] + ']}, '
            # Removing the extra comma
            txt = txt[:-2]
            txt += ']\n'

        angsposfile = '{}/atomic_angposfile.txt'.format(tetris_path)
        if os.path.exists(angsposfile):
            os.remove(angsposfile)
        with open(angsposfile, 'w') as f:
            f.write(txt)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for generating 'tetrises' (dense volumetric distributions) from volumes input.",
                                     epilog = "Example: python %(prog)s --volumes volumes --templates templates --distractors distractors --output tetrises"
                                              " --frequencies_csv Frequencies.csv --tetris_sampling_rate 16 --number_of_tetrises 48"
                                              " --dimensions 128 128 64 --density_ratio 3 --iterations 5 --insertion_distances -1 0"
                                              " --sigma 1.65 --gray_level_threshold 100 --threads 5 --grind False")

    parser.add_argument('--volumes', type=str, default='volumes', help='Default: %(default)s. Directory where volume versions of the templates and distractors exist in corresponding subdirectories.')
    parser.add_argument('--templates', type=str, default='templates', help='Default: %(default)s. Directory where all the templates PDBs exist.')
    parser.add_argument('--distractors', type=str, default='distractors', help='Default: %(default)s. Directory where all the distractors PDBs exist.')
    parser.add_argument('--output', type=str, default='volumes', help='Default: %(default)s. Directory where all the output tetrises will be stored.')
    parser.add_argument('--frequencies_csv', type=str, default='Frequenceis.csv', help='Default: %(default)s. CSV file containing all the distractors names with how many times they appear in each iteration (small distractors appear more than big ones).')

    parser.add_argument('--tetris_sampling_rate', type=int, default=16, help='Default: %(default)s. Sampling rate, must match the input volumes sampling rate.')
    parser.add_argument('--number_of_tetrises', type=int, default=48, help='Default: %(default)s. Number of tetrises to generate.')
    parser.add_argument('--dimensions', nargs=3, type=str, default=['128', '128', '64'], help=" Default: %(default)s. Dimensions for the tetris. Pass as three separate values.")
    parser.add_argument('--density_ratio', type=int, default=2, help="Default: %(default)s. Number of distractors for each template. If the template is small (say 200 kDa), use 1 or 2."
                                                                     " If the template is big (say a few MDa), use around 5.")
    parser.add_argument('--iterations', type=int, default=5, help='Default: %(default)s. Number of iterations to populate each tetris.')
    parser.add_argument('--insertion_distances', nargs=2, type=int, default=[-1, 0], help="Default: %(default)s. Using [-2, 0] generates extremely dense outputs,"
                                                                                          " if you want to reduce the density, use [0, 1] for example, or higher.")
    parser.add_argument('--sigma', type=float, default=1.65, help='Default: %(default)s. Sigma value for the lowpass filtering (for experts).')
    parser.add_argument('--gray_level_threshold', type=int, default=100, help='Default: %(default)s. Gray level threshold (for experts).')
    parser.add_argument('--threads', type=int, default=None, help='Default: %(default)s. Number of threads to use, None for all CPUs -1.')
    parser.add_argument('--grind', type=str2bool, default=False, help='Default: %(default)s. Keep trying to add molecules even when saturation is achieved.'
                                                                      ' Usually, this is not needed nor efficient, but if you want to keep trying to make'
                                                                      ' the output volumes denser, you can set it to true.')
    args = parser.parse_args()

    main(args.volumes, args.templates, args.distractors, args.output, args.frequencies_csv, args.tetris_sampling_rate,
         args.number_of_tetrises, args.dimensions, args.density_ratio, args.iterations, args.insertion_distances, args.sigma,
         args.gray_level_threshold, args.threads, args.grind)
