import numpy as np
import gc
from mpi4py import MPI
import datetime

from nbodykit.cosmology import Cosmology
from nbodykit.source.catalog import BigFileCatalog
from nbodykit.lab import FOF

# start = datetime.datetime.now()

# wdir = "../catalog/calibration/jiajun/a_0.7692/"

# dm = BigFileCatalog(wdir, dataset='1/', header='Header')
# np.savetxt("test.txt", np.array(dm['Position']), fmt="%.3f %.3f %.3f")

# end = datetime.datetime.now()
# print(end-start)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

MYPlanck = Cosmology(m_ncdm=[],
               Omega0_b=0.0491,
               Omega0_cdm=0.3156 - 0.0491, 
               h=0.6727)\
            .match(sigma8=0.831)

if rank == 0:
    start = datetime.datetime.now()

comm.barrier()

wdir = "./Nbody/"

dm = BigFileCatalog(wdir, dataset='1/', header='Header')
fof = FOF(dm, linking_length=0.200, nmin=20,)

halos = fof.to_halos(particle_mass=dm.attrs['MassTable'][1]*1e10, cosmo=MYPlanck, redshift=0.3, peakcolumn=None)
halos.save("./Nbody/", dataset="LL-0.200/", header="Header")