import argparse
import glob
import os
import random
import tqdm
import numpy as np
import prody


def main(deformation_range, deformations_per_template, templates_dir):
    templates = list(glob.glob('{}/*.pdb'.format(templates_dir)))
    nma_paths = '{}/NMA/'.format(templates_dir)

    print('Will be deforming these templates: ', templates)
    for template in tqdm.tqdm(templates):
        name = os.path.basename(template)[:-4]
        nma_vectors_path = glob.glob(nma_paths + name + '/modes/*')
        modes = []
        for vector in nma_vectors_path:
            # ignore the first six modes that correspond to rigid-body deformations and read the other modes
            if vector.endswith(('.1', '.2', '.3', '.4', '.5', '.6')):
                continue
            else:
                mode = np.loadtxt(vector)
            modes.append(mode)

        print('\nDeforming versions from ', template)
        for i in tqdm.tqdm(range(deformations_per_template)):
            pdb = prody.parsePDB(template)
            pdb_ori = prody.parsePDB(template)  # This seems unused, consider removing it if not needed
            for mode in modes:
                mode = np.array(mode)
                amplitude = (random.random() - 0.5) * 2 * deformation_range
                mode *= amplitude
                pdb._coords += mode
            new_name = '{}/{}_{}.pdb'.format(templates_dir, name, i)
            prody.writePDB(new_name, pdb)
    print('If you want to remove the deformed templates and repeat this script with different parameters,'
          ' use this command "rm templates/*centered_*.pdb"')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate deformations for input templates.',
                                     epilog="Example: python %(prog)s --working_directory templates "
                                            "--deformation_range 100 --deformations_per_template 25")

    parser.add_argument('--working_directory', type=str, default='templates', help='Default: %(default)s. The directory'
                                                                                   ' where you have the centered'
                                                                                   ' templates and their normal modes.')
    parser.add_argument('--deformation_range', type=int, default=100, help='Default: %(default)s. Range of normal mode deformations.'
                                                                           ' Usually 100 works well. If you see the '
                                                                           'molecule very deformed after excuting this step'
                                                                           ' (verify in Chimera some structures in the '
                                                                           ' working directory),'
                                                                           ' reduced it'
                                                                           ' (e.g. to 50). If the molecules are not'
                                                                           ' deformed enough, then increase it'
                                                                           ' (e.g., to 200).')
    parser.add_argument('--deformations_per_template', type=int, default=25, help='Default: %(default)s. Number of'
                                                                                  ' deformations per template. The more'
                                                                                  ' usually the better, but the more time'
                                                                                  ' needed in the next step'
                                                                                  ' (when converting to volumes).'
                                                                                  ' In the Template Learning paper, '
                                                                                  ' we tested between 25 per template'
                                                                                  ' to 100 per template. If there are '
                                                                                  ' a few templates (like 1 or 2 per '
                                                                                  ' molecule) use a higher number to '
                                                                                  ' compensate (e.g. 100 to 1000)')
    args = parser.parse_args()

    main(args.deformation_range, args.deformations_per_template, args.working_directory)