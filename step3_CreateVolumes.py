import os
import glob
import tqdm

############# USER #############
box_size = 22
voxel_size = 16
resolution = 32

############# CODE #############
os.mkdir('volumes')
os.mkdir('volumes/distractors')
os.mkdir('volumes/templates')

# converting distractors to volumes:
distractors = list(glob.glob('distractors/*.pdb'))

print('Converting distractors')
for distractor in tqdm.tqdm(distractors):
    volume_name = 'volumes/distractors/' + os.path.basename(distractor)[:-3] + 'mrc'
    args = 'e2pdb2mrc.py --box {} --res {} --apix {} {} {}'.format(box_size, resolution, voxel_size, distractor, volume_name)
    os.system(args)

templates = list(glob.glob('templates/*.pdb'))
print('Converting templates')
for template in tqdm.tqdm(templates):
    volume_name = 'volumes/templates/' + os.path.basename(template)[:-3] + 'mrc'
    args = 'e2pdb2mrc.py --box {} --res {} --apix {} {} {}'.format(box_size, resolution, voxel_size, template, volume_name)
    os.system(args)
