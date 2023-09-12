# remember to activate eman2 environment before running this script

import os
import glob
import tqdm
import yaml
import numpy


# thickness
z = 128
sirt_iterations = 0  # set to zero if you want only WBP reconstructions


def createTltfile(name):
    presorted = '{}/presorted.tlt'.format(name)
    # read the config and get starting and step
    config_file = '{}/config.yaml'.format(name)
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    start_angle = int(config['scan']['start_angle'])
    step_angle = int(config['scan']['step_angle'])
    # create and angles list
    angles = []
    for i in range(0, -start_angle + step_angle, step_angle):
        if i == 0:
            angles.append(i)
        else:
            angles.append(-i)
            angles.append(i)
    # save as a tlt file
    if os.path.exists(presorted):
        os.remove(presorted)
    numpy.savetxt(presorted, angles, fmt='%f')
    return presorted


def sortTiltSeries(tsfn, tltfn, output):
    # Unstacking the tilt series
    args = 'e2spt_tiltstacker.py --unstack {} --tltfile {}'.format(tsfn, tltfn)
    os.system(args)
    tilts = numpy.loadtxt(tltfn)
    framesfn = sorted(list(glob.glob('sptstacker_01/*.mrcs')))
    indices = sorted(range(len(tilts)), key=lambda index: tilts[index])
    sorted_tltfile = [list(tilts)[index] for index in indices]
    sorted_fns = [framesfn[index] for index in indices]

    # now creating an imod file:
    numpy.savetxt('sorted.tlt', numpy.array(sorted_tltfile), fmt="%.2f")
    with open('imod_file.txt', 'x') as f:
        f.write('{}\n'.format(len(tilts)))
        for fn in sorted_fns:
            f.write(fn + '\n')
            f.write('0\n')

    # removing the old link to the unsorted tilt series and creating the sorted tilt series stack
    command = 'newstack -fileinlist {} -output {}'.format('imod_file.txt', output)
    os.system(command)


simulation_list = glob.glob('vpp/*')
parent_dir = os.path.abspath(os.curdir)

print(simulation_list)

# sort tilt series
for simulation in tqdm.tqdm(simulation_list):
    presorted = createTltfile(simulation)
    os.chdir(simulation)
    tsfn = 'tiltseries.mrcs'
    output = 'sorted.mrcs'
    # if the file exists, delete it
    if os.path.exists('imod_file.txt'):
        os.system('rm -rf {} {}'.format('imod_file.txt', 'sptstacker*'))
    sortTiltSeries(tsfn, os.path.basename(presorted), output)
    os.chdir(parent_dir)

# # reconstruct
for simulation in tqdm.tqdm(simulation_list):
    os.chdir(simulation)
    tsfn = 'sorted.mrcs'
    reconstruction = 'sorted.rec'
    tomogram = 'tomogram.mrc'
    tltfile = 'sorted.tlt'
    args = 'tilt -InputProjections {} -OutputFile {} -TILTFILE {} -THICKNESS {}' \
           ' -FalloffIsTrueSigma 1 -RADIAL 0.5,0.5 -SHIFT 0.0,0.0 -OFFSET 0.0,0.0 -MODE 1 -PERPENDICULAR' \
           ' -AdjustOrigin -UseGPU 0 -ActionIfGPUFails 2,2'.format(tsfn, reconstruction, tltfile, z)
    if sirt_iterations is not 0:
        args += ' -FakeSIRTiterations {}'.format(sirt_iterations)
    os.system(args)
    args = 'trimvol -rx {} {}'.format(reconstruction, tomogram)
    os.system(args)
    os.remove(reconstruction)
    os.chdir(parent_dir)

