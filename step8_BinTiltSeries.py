import os
import glob
import tqdm

# binning pixel size
PixelSize = 8  # in Angstroms, has to be integer for imod to accept it

simulation_list = glob.glob('parakeet/*')
parent_dir = os.path.abspath(os.curdir)

print(simulation_list)

for simulation in tqdm.tqdm(simulation_list):
    os.chdir(simulation)
    tilt_series_h5 = 'image.h5'
    tilt_series_mrc = 'image.mrc'
    if os.path.exists(tilt_series_mrc):
        os.remove(tilt_series_mrc)
    # export the data to mrc:
    os.system('parakeet.export {} -o {}'.format(tilt_series_h5, tilt_series_mrc))
    # bin the data with imod
    tilt_series_binned = 'tiltseries.mrcs'
    if os.path.exists(tilt_series_binned):
        os.remove(tilt_series_binned)
    os.system('newstack -input {} -output {} -taper 1,1 -bin {} '
              '-antialias 6 -imagebinned 1.0'.format(tilt_series_mrc, tilt_series_binned, PixelSize))
    # remove the original size tilt series
    os.system('rm {}'.format(tilt_series_mrc))
    os.chdir(parent_dir)
