import mrcfile
import glob

tomolist = glob.glob('*.mrc')

for tomo in tomolist:
    with mrcfile.open(tomo) as mrc:
        newdata = -mrc.data
    mrcfile.write(tomo, newdata, overwrite=True)
