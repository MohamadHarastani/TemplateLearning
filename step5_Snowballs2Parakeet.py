import os
import glob
import tqdm

from scipy.spatial.transform import Rotation

snowball_sampling_rate = 16

snowballs_path = glob.glob('snowballs/*')
templates = glob.glob('templates/*')
distractors = glob.glob('distractors/*')

templates_base_names = [os.path.basename(template)[:-4] for template in templates]
distractors_base_names = [os.path.basename(distractor)[:-4] for distractor in distractors]


for snowball_path in tqdm.tqdm(snowballs_path):
    coordinates_tables = glob.glob('{}/coordinates/*.txt'.format(snowball_path))
    angles_tables = glob.glob('{}/angles/*.txt'.format(snowball_path))
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
            coords = [float(i) * snowball_sampling_rate for i in coords.split(',')]
            txt += '{position: [' + '{}'.format(coords)[1:-1]
            # converting the angles to rotation vectors
            angs = [float(i) for i in angs.split(',')]
            T = Rotation.from_euler('ZYZ', [angs[0], angs[1], angs[2]], degrees=True).inv().as_matrix()
            angs = list(Rotation.from_matrix(T).as_rotvec())
            txt += '], orientation: [' + '{}'.format(angs)[1:-1] + ']}, '
        # Removing the extra comma
        txt = txt[:-2]
        txt += ']\n'

    angsposfile = '{}/atomic_angposfile.txt'.format(snowball_path)
    if os.path.exists(angsposfile):
        os.remove(angsposfile)
    with open(angsposfile, 'w') as f:
        f.write(txt)
