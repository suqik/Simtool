import sys
import numpy as np
import pyccl as ccl
from decimal import Decimal, ROUND_HALF_UP

use_message = "Usage: python make-redshifts.py [options] [value]. \n\
Options: \n\
    -h --help        This message\n\
    -d --dist        single distance in Mpc/h\n\
    -l --dist_list   distance list in Mpc/h\n\
                     example: \"[10,20,30,40]\"\n\
    -a --dist_array  linear distance array (start, end, step) in Mpc/h\n\
                     example: \"(0 1000 10)\"\n\
    -hubble          hubble parameter in km/s/Mpc. Default is 0.7\n\
    -Om              Total matter density parameter. Default is 0.3, or 1-OL\n\
    -OL              Dark Energy density parameter. Default is 0.7, or 1-Om"

if len(sys.argv) < 2 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
    print(use_message)
    exit(1)

h  = 0.7
Om = None
OL = None
singlef = 0
dlistf = 0
darrf  = 0
for idx, param in enumerate(sys.argv):
    if param == '-hubble':
        h = float(sys.argv[idx+1])
    if param == '-Om':
        Om = float(sys.argv[idx+1])
    if param == '-OL':
        OL = float(sys.argv[idx+1])
    if param == '-d' or param == "--dist":
        singlef = 1
        d = float(sys.argv[idx+1])
    if param == '-l' or param == "--dist_list":
        dlistf = 1
        dlist = str(sys.argv[idx+1])[1:-1].split(",")
    if param == '-a' or param == "--dist_array":
        darrf = 1
        start, end, step = str(sys.argv[idx+1])[1:-1].split(" ")
        darr = np.arange(float(start), float(end), float(step))

if (singlef+dlistf+darrf) != 1:
    print("Only 1 type of distance should be given.")
    print(use_message)
    exit(1)

if not Om and not OL:
    Om = 0.3
    OL = 0.7
elif not Om:
    Om = 1 - OL
elif not OL:
    OL = 1 - Om

cosmo = ccl.Cosmology(h=h, Omega_b=0.05, Omega_c=Om-0.05, n_s=1.0, sigma8=0.831)

if singlef:
    z = 1./(ccl.scale_factor_of_chi(cosmo, d/h))-1
    print('{'+"{:.8f}".format(z)+'}')
if dlistf:
    zlist = []
    for di in dlist:
        zi = 1./ccl.scale_factor_of_chi(cosmo, float(di)/h)-1
        zlist.append(str(Decimal(str(zi)).quantize(Decimal("0.00000001"), ROUND_HALF_UP)))
    s = '{' + ', '.join(zlist) + '}'
    print(s)
if darrf:
    zlist = []
    for i in range(len(darr)):
        zi = 1./ccl.scale_factor_of_chi(cosmo, float(darr[i])/h)-1
        zlist.append(str(Decimal(str(zi)).quantize(Decimal("0.00000001"), ROUND_HALF_UP)))
    s = '{' + ', '.join(zlist) + '}'
    print(s)