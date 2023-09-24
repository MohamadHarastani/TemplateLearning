import glob
import os
import prody
import tqdm

# original templates, centered, with all atoms
templates_original_centered = "templates"

# coarse grained templates with their NMA done in ProDy
templates_cg = "templates_CG_NMA"

template_list = glob.glob('{}/*.pdb'.format(templates_original_centered))

for template in tqdm.tqdm(template_list):
    template = os.path.basename(template[:-4])

    output_path = '{}/NMA/{}/modes'.format(templates_original_centered, template)
    os.makedirs(output_path)

    full = prody.parsePDB('{}/{}.pdb'.format(templates_original_centered, template))

    coarse_grained = full.select('name CA P')
    coarse_modes = prody.parseNMD('{}/NMA/{}/modes.nmd'.format(templates_cg, template))
    id = 1
    for mode in coarse_modes[0]:
        id_path = '{}/{}'.format(output_path, id)
        os.makedirs(id_path)
        newmode = prody.extendMode(mode, coarse_grained, full, norm=True)
        prody.writeScipionModes(id_path, newmode[0])
        if id == 1:
            prody.writePDB('{}/{}_new.pdb'.format(templates_original_centered, template), newmode[1])
        os.system('mv {}/modes/vec.1 {}/vec.{}'.format(id_path, output_path, id))
        os.system('rm -rf {}'.format(id_path))
        id += 1
