# remember to activate your eman2 conda environment before launching this script
import os
import glob
import tqdm
import numpy
import random

# paramters to set
number_of_snowballs = 100
dimentions = ['256', '256', '64']
iterations = 5
insersion_distances = [1, 2]
sigma = 2
gray_level_threshold = 100
threads = None
#

templates = list(glob.glob('volumes/templates/*.mrc'))
distractors_and_frequencies = numpy.loadtxt('Frequencies.csv', delimiter=',', dtype='str')
distractors_and_frequencies = distractors_and_frequencies.transpose()
os.mkdir('snowballs/')

for snowball_number in tqdm.tqdm(range(number_of_snowballs)):
    molecules_list = []
    snowball_path = 'snowballs/{}'.format(snowball_number)
    coordinates_tables_path = '{}/coordinates'.format(snowball_path)
    angles_tables_path = '{}/angles'.format(snowball_path)
    output_volume_path = '{}/output_volume.mrc'.format(snowball_path)
    os.mkdir(snowball_path)
    os.mkdir(coordinates_tables_path)
    os.mkdir(angles_tables_path)

    random.shuffle(templates)
    # fix the names to include the path and _centered.mrc, and add the templates
    for i in range(len(distractors_and_frequencies)):
        new_line = ['volumes/distractors/{}_centered.mrc'.format(distractors_and_frequencies[i, 0]), distractors_and_frequencies[i, 1]]
        molecules_list.append(new_line)
        new_line = [templates[i], distractors_and_frequencies[i, 1]]
        molecules_list.append(new_line)
    # now run the snowball
    molecules = list(numpy.array(molecules_list)[:,0])
    frequencies = list(numpy.array(molecules_list)[:,1])
    arg = "python functions/multisnowball_python.py --mol '{}' --dim '{}' --frequencies '{}' --iterations {}" \
          " --coordinate_table {} --angle_table {} --output_volume {} --insersion_distances '{}' --sigma {}" \
          " --gray_level_threshold {} --threads {}".format(molecules, dimentions, frequencies, iterations,
                                                           coordinates_tables_path, angles_tables_path,
                                                           output_volume_path, insersion_distances, sigma,
                                                           gray_level_threshold, threads)
    os.system(arg)
    # print(arg)