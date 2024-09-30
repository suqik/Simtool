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

def mk_mesh_cat(pos, Nmesh, boxsize, MAS="NGP"):
    delta = np.zeros((Nmesh, Nmesh, Nmesh), dtype=np.float32)
    pos = pos.astype(np.float32)

    MASL.MA(pos, delta, boxsize, MAS)
    del pos
    dr = boxsize/float(Nmesh)
    mesh_pos = np.linspace(dr/2., boxsize-dr/2., Nmesh)
    mesh_x, mesh_y, mesh_z = np.meshgrid(mesh_pos, mesh_pos, mesh_pos, indexing="ij")

    grid_non_zero = (delta != 0)
    mesh_cat_pos = np.c_[mesh_x[grid_non_zero], mesh_y[grid_non_zero], mesh_z[grid_non_zero]]
    mesh_cat_weight = delta[grid_non_zero]

    return mesh_cat_pos, mesh_cat_weight

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