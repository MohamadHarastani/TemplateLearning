import argparse
import yaml
import os
import numpy as np
import glob

# a function to accept boolean values as input to the parser
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# a function to create a tlt file using the yaml config file of parakeet
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
    np.savetxt(presorted, angles, fmt='%f')
    return presorted

# a function that sorts a tilt series
def sortTiltSeries(tsfn, tltfn, output):
    # Unstacking the tilt series
    args = 'e2spt_tiltstacker.py --unstack {} --tltfile {}'.format(tsfn, tltfn)
    os.system(args)
    tilts = np.loadtxt(tltfn)
    framesfn = sorted(list(glob.glob('sptstacker_01/*.mrcs')))
    indices = sorted(range(len(tilts)), key=lambda index: tilts[index])
    sorted_tltfile = [list(tilts)[index] for index in indices]
    sorted_fns = [framesfn[index] for index in indices]

    # now creating an imod file:
    np.savetxt('sorted.tlt', np.array(sorted_tltfile), fmt="%.2f")
    with open('imod_file.txt', 'x') as f:
        f.write('{}\n'.format(len(tilts)))
        for fn in sorted_fns:
            f.write(fn + '\n')
            f.write('0\n')

    # removing the old link to the unsorted tilt series and creating the sorted tilt series stack
    command = 'newstack -fileinlist {} -output {}'.format('imod_file.txt', output)
    os.system(command)


# a function to resize 3D to new dimensions
def resize(data, x, y, z):
    """Resize the data.

    @param data: input data.
    @param x: resized dimension x.
    @param y: resized dimension y.
    @param z: resized dimension z.

    @return: resized data.
    """
    s = data.shape
    from scipy import mgrid, array
    from scipy.ndimage import map_coordinates
    grid = mgrid[0:s[0] - 1:x * 1j, 0:s[1] - 1:y * 1j, 0:s[2] - 1:z * 1j]
    d = map_coordinates(data, grid, order=0)
    return d