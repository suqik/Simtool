import sys
import numpy as np
from nbodykit.source.catalog.file import BigFileCatalog

file = sys.argv[1]
ftype = sys.argv[2]
ofile = sys.argv[3]

if ftype.lower() == 'dm':
    dataset = '1/'
elif ftype.lower() == 'halo':
    dataset = 'RFOF/'
else:
    print("Cannot recognize the dataset type. Only `dm` or `halo` (Case insensitive).")
    exit(1)

cat = BigFileCatalog(file, dataset=dataset, header='Header')
np.savetxt(ofile, np.array(cat['Position']), fmt='%.3f %.3f %.3f')