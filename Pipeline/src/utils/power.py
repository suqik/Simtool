import numpy as np
import MAS_library as MASL
import Pk_library as PKL

def cat2mesh(pos, Nmesh, boxsize, MAS="CIC"):
    delta = np.zeros((Nmesh, Nmesh, Nmesh), dtype=np.float32)
    pos = pos.astype(np.float32)

    MASL.MA(pos, delta, boxsize, MAS)
    delta /= np.mean(delta, dtype=np.float64)
    delta -= 1.0

    return delta

def get_powerspec_cat(pos, Nmesh, boxsize, MAS="CIC", kmin=None, kmax=None, nk=None, multipoles=[0], threads=None):
    delta = cat2mesh(pos, Nmesh, boxsize, MAS)
    Pk = PKL.Pk(delta, boxsize, axis=0, MAS=MAS, threads=threads)

    k = Pk.k3D

    if kmin is None or kmin < k[0]:
        kmin = k[0]
    if kmax is None or kmax > k[-1]:
        kmax = k[-1]
    if nk is None:
        nk = len(k)

    kout = np.logspace(np.log10(kmin), np.log10(kmax), nk)
    Pkout = []
    
    for ell in multipoles:
        Pk_interp = np.interp(np.log10(kout), np.log10(k), Pk.Pk[:,int(ell//2)])
        Pkout.append(Pk_interp)

    return kout, Pkout