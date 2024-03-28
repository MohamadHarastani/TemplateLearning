import argparse
import glob
import os
import multiprocessing
import tqdm
import prody
import numpy as np
from functions.utilities import str2bool


def main(input, output, center_templates, coarse_grain, course_grain_command, cif, modes_number, cutoff,
         spring_constant, threads, debug):
    # create the output directory if it does not exist:
    if not os.path.exists(output):
        os.mkdir(output)

    # in is a temporary path to handle coarse graining
    backup_path = '{}/backups'.format(output)
    if not os.path.exists(backup_path):
        os.mkdir(backup_path)

    if cif:
        templates = list(glob.glob('{}/*.cif'.format(input)))
    else:
        templates = list(glob.glob('{}/*.pdb'.format(input)))

    if center_templates:
        print('Centering the templates')
        for template in tqdm.tqdm(templates):
            pdb = prody.parsePDB(template)
            pdb._coords -= np.mean(pdb._coords, axis=1)
            new_name = '{}/{}'.format(output, os.path.basename(template)[:-4] + '_centered')
            prody.writePDB(new_name, pdb)

    # if the templates won't get centered, then we still need to make sure the templates are readible by Prody
    # and save them as PDB files if they were originally cif
    else:
        for template in tqdm.tqdm(templates):
            if cif:
                templates = list(glob.glob('{}/*.cif'.format(input)))
            else:
                templates = list(glob.glob('{}/*.pdb'.format(input)))
            pdb = prody.parsePDB(template)
            new_name = '{}/{}'.format(output, os.path.basename(template)[:-4])
            prody.writePDB(new_name, pdb)

    # Use the centered copies or just the PDB verified copies of the templates in the rest of the protocol
    templates = list(glob.glob('{}/*.pdb'.format(output)))

    # The calculated normal modes will be placed here:
    nma_paths = '{}/NMA'.format(output)
    os.mkdir(nma_paths)

    print('Calculating NMA for each template')
    for template in tqdm.tqdm(templates):
        name = os.path.basename(template)[:-4]
        nma_path = nma_paths + '/' + name
        os.mkdir(nma_path)
        if coarse_grain:
            backup_template_name = template
            cg_template = '{}_cg.pdb'.format(template[:-4])
            command = "prody select '{}' {} -o {}".format(course_grain_command, template, cg_template)
            os.system(command)
            template = cg_template

        # for details, see Prody documentation for ANM NMA
        args = 'prody anm {0} -s "all" --altloc "all"  --hessian --export-scipion --npzmatrices ' \
               '--npz -o {1} -p modes -n {2} -g {3} -c "{4}" -P {5} --turbo'.format(template,
                                                                                    nma_path,
                                                                                    modes_number,
                                                                                    spring_constant,
                                                                                    cutoff,
                                                                                    threads)
        os.system(args)

        if coarse_grain:
            # rename previous normal mode path to cg_
            os.system('mv {} {}_cg'.format(nma_path, nma_path))
            # make a new normal mode path
            os.mkdir(nma_path)
            os.mkdir('{}/modes'.format(nma_path))

            # coarse grained templates will also be moved to the backup
            full = prody.parsePDB(backup_template_name)
            coarse_grained = full.select(course_grain_command)
            coarse_modes = prody.parseNMD('{}_cg/modes.nmd'.format(nma_path))
            id = 1
            for mode in coarse_modes[0]:
                id_path = '{}/{}'.format(nma_path, id)
                os.makedirs(id_path)
                # extending the modes (interpolating) to the full atomic strucutre
                newmode = prody.extendMode(mode, coarse_grained, full, norm=True)
                prody.writeScipionModes(id_path, newmode[0])
                if id == 1:
                    new_template_name = '{}_new.pdb'.format(template[:-4])
                    prody.writePDB(new_template_name, newmode[1])
                # Since each mode has to be saved separately after interpolating (as vec.1), so adapting the numbering
                os.system('mv {}/modes/vec.1 {}/modes/vec.{}'.format(id_path, nma_path, id))
                if not debug:
                    os.system('rm -rf {}'.format(id_path))
                id += 1

            if not debug:
                os.system('rm -rf {}_cg'.format(nma_path))

            os.system('mv {} {}'.format(backup_template_name, backup_path))
            os.system('mv {} {}'.format(template, backup_path))
            os.system('mv {} {}'.format(new_template_name, backup_template_name))

    if not debug:
        os.system('rm -rf {}'.format(backup_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='This script takes input atomic structures placed in the same directory'
                    ' and does the following :'
                    '1- Centers them to the origin of coordinates (important! but optional if its already done). '
                    '2- It will (optionally) coarse grain them to CA and P atoms (important to large structures for speed). '
                    '3- It will perform Normal Mode Analysis to each template'
                    ' (or coarse grained version of the template). '
                    '4- In case they are coarse grained, it will interpolate the output'
                    ' modes to the full structure. '
                    '5- The output of this script is: '
                    'a- centered full atomic structures. '
                    'b- corresponding normal modes.',
        epilog="Example: python %(prog)s --input_directory input_templates --output_directory templates"
               " --center_templates True --cif False --coarse_grain True --course_grain_command name CA P "
               " --modes_number 20 --cutoff 15 --spring_constant 1 --keep_intermediate_files False")


    parser.add_argument('--input_directory', type=str, default='input_templates',
                        help='Default: %(default)s. Directory for input atomic structures to be used as templates.')
    parser.add_argument('--output_directory', type=str, default='templates',
                        help='Default: %(default)s. Path to save pre-processed templates and normal modes.')
    parser.add_argument('--center_templates', type=str2bool, default=True,
                        help='Default: %(default)s. Center templates before processing')
    parser.add_argument('--cif', type=str2bool, default=False,
                        help='Templates are in cif instead of pdb format (use CIF either or PDB and not a mix). Default: %(default)s')
    parser.add_argument('--coarse_grain', type=str2bool, default=True,
                        help='Apply coarse graining to templates. Default: %(default)s')
    parser.add_argument('--course_grain_command', type=str, default="name CA P",
                        help='ProDy selection command for coarse graining (will be ignored if --coarse_grain is False).'
                             ' See Prody documentation if you want other coarse graining action. Default: %(default)s')
    parser.add_argument('--modes_number', type=int, default=20,
                        help='Number of modes to calculate. Default: %(default)s')
    parser.add_argument('--cutoff', type=int, default=15,
                        help='Cutoff distance for ANM NMA calculations. Default: %(default)s')
    parser.add_argument('--spring_constant', type=int, default=1,
                        help='Spring constant for ANM calculations. Default: %(default)s')
    parser.add_argument('--threads', type=int, default=multiprocessing.cpu_count() - 1,
                        help='Number of threads to use. Default: %(default)s')
    parser.add_argument('--keep_intermediate_files', type=str2bool, default=True,
                        help='Keep coarse grained modes on the disk for inspection. Default: %(default)s')

    args = parser.parse_args()

    main(args.input_directory, args.output_directory, args.center_templates, args.coarse_grain, args.course_grain_command, args.cif,
         args.modes_number, args.cutoff, args.spring_constant, args.threads, args.keep_intermediate_files)
