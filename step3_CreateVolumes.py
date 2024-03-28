import argparse
import os
import glob
import tqdm


def main(templates, distractors, output, box_size, voxel_size, resolution):
    if not os.path.exists(output):
        os.makedirs(output)
    if not os.path.exists('{}/distractors'.format(output)):
        os.makedirs('{}/distractors'.format(output))
    if not os.path.exists('{}/templates'.format(output)):
        os.makedirs('{}/templates'.format(output))
        
        
    # Converting templates to volumes:
    templates = list(glob.glob('{}/*.pdb'.format(templates)))
    print('Converting templates\n')
    for template in tqdm.tqdm(templates):
        volume_name = '{}/templates/'.format(output) + os.path.basename(template)[:-3] + 'mrc'
        args = 'e2pdb2mrc.py --box {} --res {} --apix {} {} {}'.format(box_size, resolution, voxel_size, template, volume_name)
        os.system(args)

    # Converting distractors to volumes:
    distractors = list(glob.glob('{}/*.pdb'.format(distractors)))
    print('Converting distractors\n')
    for distractor in tqdm.tqdm(distractors):
        volume_name = '{}/distractors/'.format(output) + os.path.basename(distractor)[:-3] + 'mrc'
        args = 'e2pdb2mrc.py --box {} --res {} --apix {} {} {}'.format(box_size, resolution, voxel_size, distractor, volume_name)
        os.system(args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert PDB files to MRC volumes with specified parameters.',
                                     epilog="Example: python %(prog)s --templates templates --distractors distractors"
                                            " --output volumes --box_size 22 --voxel_size 16 --resolution 32")
    parser.add_argument('--templates', type=str, default='templates', help='Default: %(default)s. Directory where all the templates PDBs exist.')
    parser.add_argument('--distractors', type=str, default='distractors', help='Default: %(default)s. Directory where all the distractors PDBs exist.')        
    parser.add_argument('--output', type=str, default='volumes', help='Default: %(default)s. Directory where all the output volumes will be stored.')
    parser.add_argument('--box_size', type=int, default=22, help='Box size for the conversion. Default: %(default)s.')
    parser.add_argument('--voxel_size', type=int, default=16, help='Voxel size for the conversion. Default: %(default)s.')
    parser.add_argument('--resolution', type=int, default=32, help='Resolution for the conversion. Default: %(default)s.')

    args = parser.parse_args()

    main(args.templates, args.distractors, args.output, args.box_size, args.voxel_size, args.resolution)
    
