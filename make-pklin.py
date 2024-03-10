import sys
import numpy

ofile = "powerspec.txt"

if sys.argv[1] == '-h' or sys.argv[1] == '--help':
      print("Usage: python make-pklin.py [options]\n\
Options:\n\
   -h --help     This message.\n\
   -hubble       hubble parameter in km/s/Mpc\n\
   -OmegaM       total matter density parameter\n\
   -Omegab       baryon density parameter\n\
   -sigma8       matter fluctuation at 8Mpc/h scale\n\
   -ofile        path of output file.")
      exit(1)

for idx, param in enumerate(sys.argv[1:]):
   if param[1:] == 'hubble':
      h = float(sys.argv[idx+2])
   elif param[1:] == 'OmegaM':
      Omega0_m = float(sys.argv[idx+2])
   elif param[1:] == 'Omegab':
      Omega0_b = float(sys.argv[idx+2])
   elif param[1:] == 'sigma8':
      sigma8 = float(sys.argv[idx+2])
   elif param[1:] == 'ofile':
      ofile = (sys.argv[idx+2])

from nbodykit.cosmology import LinearPower, Cosmology
from nbodykit.lab import *

MYPlanck = Cosmology(m_ncdm=[],
               Omega0_b=Omega0_b,
               Omega0_cdm=Omega0_m - Omega0_b, 
               h=h)\
            .match(sigma8=sigma8)

pklin0 = LinearPower(MYPlanck, redshift=0)

k = numpy.logspace(-3, 2, 10000, endpoint=True)

numpy.savetxt(ofile, list(zip(k, pklin0(k))))

