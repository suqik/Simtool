import sys
import numpy

ofile = "powerspec.txt"

if len(sys.argv) != 13:
      print("Usage: python make-pklin.py [options]\n\
Options:\n\
   -h --help     This message.\n\
   -z            redshift of power spectrum\n\
   -hubble       hubble parameter in km/s/Mpc\n\
   -OmegaM       total matter density parameter\n\
   -Omegab       baryon density parameter\n\
   -sigma8       matter fluctuation at 8Mpc/h scale\n\
   -ofile        path of output file.")
      exit(1)

for idx, param in enumerate(sys.argv):
   if param[1:] == 'z':
      if float(sys.argv[idx+1]) < 0:
         raise ValueError("Input redshift must larger than 0!")
      else:
         z = float(sys.argv[idx+1])
   if param[1:] == 'hubble':
      h = float(sys.argv[idx+1])
   elif param[1:] == 'OmegaM':
      Omega0_m = float(sys.argv[idx+1])
   elif param[1:] == 'Omegab':
      Omega0_b = float(sys.argv[idx+1])
   elif param[1:] == 'sigma8':
      sigma8 = float(sys.argv[idx+1])
   elif param[1:] == 'ofile':
      ofile = (sys.argv[idx+1])

from nbodykit.cosmology import LinearPower, Cosmology
from nbodykit.lab import *

MYPlanck = Cosmology(m_ncdm=[],
               Omega0_b=Omega0_b,
               Omega0_cdm=Omega0_m - Omega0_b, 
               h=h)\
            .match(sigma8=sigma8)

pklin0 = LinearPower(MYPlanck, redshift=z)

k = numpy.logspace(-3, 2, 10000, endpoint=True)

numpy.savetxt(ofile, list(zip(k, pklin0(k))))

from matplotlib import pyplot as plt
plt.loglog(k, pklin0(k))
plt.xlabel("k [h/Mpc]")
plt.ylabel("Pk [(Mpc/h)^3]")
plt.savefig("pklin-z100.png", format='png', dpi=400, bbox_inches='tight')