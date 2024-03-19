import os
import glob
import tqdm
import numpy
import random
from scipy.spatial.transform import Rotation

############# USER #############
tetris_sampling_rate = 16  # Attention: has to match the volumes sampling rate
number_of_tetrises = 48
dimentions = ['128', '128', '64']
iterations = 5
insersion_distances = [-2, 0]
sigma = 1.65
gray_level_threshold = 100
threads = None  # None is ( all CPUs -1 ) , if you want to use less CPUs, set it here
grind = False  # grind will try to add more molecules, but the tetris dimention should be big enough for it to work, keep it False unless you are an expert
density_ratio = 2  #  number of distractors for each template. Set it to 1 if your template is too small, 2 or more if the tempalte is big

############# CODE #############
templates = list(glob.glob('volumes/templates/*.mrc'))
# repleat the templates just in case a template is alone, or the number of templates are less than distractors
templates = templates * 200

distractors_and_frequencies = numpy.loadtxt('Frequencies.csv', delimiter=',', dtype='str')
distractors_and_frequencies = distractors_and_frequencies.transpose()
os.mkdir('tetrises/')

for tetris_number in tqdm.tqdm(range(number_of_tetrises)):
    molecules_list = []
    tetris_path = 'tetrises/{}'.format(tetris_number)
    coordinates_tables_path = '{}/coordinates'.format(tetris_path)
    angles_tables_path = '{}/angles'.format(tetris_path)
    output_volume_path = '{}/output_volume.mrc'.format(tetris_path)
    os.mkdir(tetris_path)
    os.mkdir(coordinates_tables_path)
    os.mkdir(angles_tables_path)

    random.shuffle(templates)
    # fix the names to include the path and _centered.mrc, and add the templates
    for i in range(len(distractors_and_frequencies)):
        new_line = ['volumes/distractors/{}_centered.mrc'.format(distractors_and_frequencies[i, 0]),
                    distractors_and_frequencies[i, 1]]
        molecules_list.append(new_line)
        if i % density_ratio == 0:
            new_line = [templates[i], distractors_and_frequencies[i, 1]]

        molecules_list.append(new_line)

    # now run the tetris
    molecules = list(numpy.array(molecules_list)[:, 0])
    frequencies = list(numpy.array(molecules_list)[:, 1])
    arg = "python functions/multitetris_python.py --mol '{}' --dim '{}' --frequencies '{}' --iterations {}" \
          " --coordinate_table {} --angle_table {} --output_volume {} --insersion_distances '{}' --sigma {}" \
          " --gray_level_threshold {} --threads {} --grind {}".format(molecules, dimentions, frequencies, iterations,
                                                           coordinates_tables_path, angles_tables_path,
                                                           output_volume_path, insersion_distances, sigma,
                                                           gray_level_threshold, threads, grind)
    os.system(arg)

# converting the tetrises coordinates to parakeet coordinates
tetrises_path = glob.glob('tetrises/*')
templates = glob.glob('templates/*')
distractors = glob.glob('distractors/*')

templates_base_names = [os.path.basename(template)[:-4] for template in templates]
distractors_base_names = [os.path.basename(distractor)[:-4] for distractor in distractors]


for tetris_path in tqdm.tqdm(tetrises_path):
    coordinates_tables = glob.glob('{}/coordinates/*.txt'.format(tetris_path))
    angles_tables = glob.glob('{}/angles/*.txt'.format(tetris_path))
    txt = ''
    for coordinates_table, angles_tables in zip(coordinates_tables, angles_tables):
        basename = os.path.basename(coordinates_table)[:-4]
        if basename in templates_base_names:
            filename = '../../templates/{}.pdb'.format(basename)
        elif basename in distractors_base_names:
            filename = '../../distractors/{}.pdb'.format(basename)
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
