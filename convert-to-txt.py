import sys
import numpy as np
import warnings
from nbodykit.source.catalog.file import BigFileCatalog

param_list = ['fclass', 'fname', 'ftype', 'cut', 'ofile']

# fclass = sys.argv[1]
# file = sys.argv[2]
# ftype = sys.argv[3]
# cut = sys.argv[4]
# ofile = sys.argv[5]

### load input params
for idx, param in enumerate(sys.argv):
    if param[0] == '-':
        param = param[1:]
        value = sys.argv[idx+1]
        if param not in param_list:
            warnings.warn(f"Cannot recognize parameter `{param}`. Will be ignored.")
        if param == 'fclass':
            fclass = value
        if param == 'fname':
            file = value
        if param == 'ftype':
            ftype = value
        if param == 'cut':
            LINF = True
            RINF = True
            minv, maxv = value.split(' ')
            if 'inf' not in minv.lower():
                if float(minv) < 0:
                    warnings.warn("Mass cannot smaller than 0. Set min M to be 0.")
                else:
                    LINF = False
                    minv = float(minv)
            if 'inf' not in maxv.lower():
                RINF = False
                maxv = float(minv)
        if param == 'ofile':
            ofile = value

if fclass.lower() == 'bigfile':
    if ftype.lower() == 'dm':
        dataset = '1/'
    elif ftype.lower() == 'halo_fof':
        dataset = 'LL-0.200/'
    elif ftype.lower() == 'halo_rfof':
        dataset = 'RFOF/'
    else:
        print("Cannot recognize the dataset type. Only `dm` or `halo_fof` or `halo_rfof` (Case insensitive).")
        exit(1)

    cat = BigFileCatalog(file, dataset=dataset, header='Header')
    if not LINF or not RINF:
        if ftype.lower() == 'halo_rfof':
            mass = cat.attrs['M0'][0]*(cat['Length']).compute()*1e10
        elif ftype.lower() == 'halo_fof':
            mass = cat['Mass'].compute()
        if not LINF and RINF:
            cut = (mass > minv)
        if LINF and not RINF:
            cut = (mass < maxv)
        if not LINF and not RINF:
            cut = ((mass > minv) & (mass < maxv))
    else:
        cut = ...
    if cat.attrs['UnitLength_in_cm'][0] - 3.08568e+21 < 1e-5:
        cat['Position'] = cat['Position']/1e3
    Nhalo = len(((cat['Position'])[cut]))
    f = open(ofile, "w+", encoding="utf-8")
    f.write(f"# Nh = {Nhalo}\n")
    np.savetxt(ofile, ((cat['Position'])[cut]).compute(), fmt='%.3f %.3f %.3f')
    f.close()

if fclass.lower() == 'ahf':
    cat = np.loadtxt(file)
    if not LINF or not RINF:
        mass = cat[:,3]
        if not LINF and RINF:
            cut = (mass > minv)
        if LINF and not RINF:
            cut = (mass < maxv)
        if not LINF and not RINF:
            cut = ((mass > minv) & (mass < maxv))
    else:   
        cut = ...
    Nhalo = len(cat[:,5:8][cut])
    f = open(ofile, "w+", encoding="utf-8")
    f.write(f"# Nh = {Nhalo}\n")
    np.savetxt(ofile, cat[:,5:8][cut]/1000., fmt='%.3f %.3f %.3f')
    f.close()

