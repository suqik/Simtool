import sys
import numpy as np
from nbodykit.source.catalog.file import BigFileCatalog

fclass = sys.argv[1]
file = sys.argv[2]
ftype = sys.argv[3]
ofile = sys.argv[4]

if fclass.lower() == 'bigfile':
    if ftype.lower() == 'dm':
        dataset = '1/'
    elif ftype.lower() == 'halo':
        dataset = 'RFOF/'
    else:
        print("Cannot recognize the dataset type. Only `dm` or `halo` (Case insensitive).")
        exit(1)

    cat = BigFileCatalog(file, dataset=dataset, header='Header')
    np.savetxt(ofile, np.array(cat['Position']), fmt='%.3f %.3f %.3f')
if fclass.lower() == 'ahf':
    cat = np.loadtxt(file)
    np.savetxt(ofile, cat[:,5:8]/1000., fmt='%.3f %.3f %.3f')
