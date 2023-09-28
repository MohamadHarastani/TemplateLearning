import glob
import os
import multiprocessing
import tqdm
import prody
import numpy as np

############# USER #############
# Templates should be placed in the folder "templates" and should be either in cif or pdb format (don't mix the two)
center_templates = True  # set True if you want your templates to be centered, otherwise False
coarse_grain = True  # set to True to go faster and coarse grain the templates before
course_grain_command = "name CA P"  # selecting only CA and P atoms while calculating normal modes. Will be ignored if course_grain is set ot False. If you want to select something else, check the ProDy documentation.
cif = False  # set True if your templates are in cif format, otherwise False
# general parameters (keep default unless you're an expert)
modes_number = 20
cutoff = 15
spring_constant = 1
threads = multiprocessing.cpu_count() - 1

############# CODE #############
backup_path = 'templates/backup_templates'
os.mkdir(backup_path)

if center_templates:
    # the centered tempalte will be saved as PDB if was in cif format
    if cif:
        templates = list(glob.glob('templates/*.cif'))
    else:
        templates = list(glob.glob('templates/*.pdb'))

    for template in tqdm.tqdm(templates):
        pdb = prody.parsePDB(template)
        pdb._coords -= np.mean(pdb._coords, axis=1)
        new_name = template[:-4] + '_centered'
        prody.writePDB(new_name, pdb)
        os.system('mv {} {}'.format(template, backup_path))

templates = list(glob.glob('templates/*.pdb'))

# This script will create random varations of your templates using normal mode analysis (NMA)
# The calculated normal modes will be placed here:
nma_paths = 'templates/NMA'
os.mkdir(nma_paths)

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
        os.system(command)

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

        # make a path where the old templates will be placed

        full = prody.parsePDB(backup_template_name)
        coarse_grained = full.select(course_grain_command)

        coarse_modes = prody.parseNMD('{}_cg/modes.nmd'.format(nma_path))
        id = 1
        for mode in coarse_modes[0]:
            id_path = '{}/{}'.format(nma_path, id)
            os.makedirs(id_path)
            newmode = prody.extendMode(mode, coarse_grained, full, norm=True)
            prody.writeScipionModes(id_path, newmode[0])
            if id == 1:
                new_template_name = '{}_new.pdb'.format(template[:-4])
                prody.writePDB(new_template_name, newmode[1])
            os.system('mv {}/modes/vec.1 {}/modes/vec.{}'.format(id_path, nma_path, id))
            os.system('rm -rf {}'.format(id_path))
            id += 1

        os.system('mv {} {}'.format(backup_template_name, backup_path))
        os.system('mv {} {}'.format(template, backup_path))
        os.system('mv {} {}'.format(new_template_name, backup_template_name))
