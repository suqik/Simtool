import os
import numpy as np
from astropy.table import Table
import heapq

rockstar_idx_dict = {
    "pos": [8,9,10],
    "mvir": 2,
    "vel": [11,12,13],
    "x": 8,
    "y": 9,
    "z": 10,
    "vx": 11,
    "vy": 12,
    "vz": 13
}

def cal_zrsd(OmegaM, redshift, r_los, vel_los):
    hubble = 100 * np.sqrt(OmegaM * (redshift + 1.0)**3 + (1-OmegaM))
    rsd_shift = (redshift + 1.0) / hubble
    s_los = r_los + vel_los*rsd_shift
    return s_los

def load_halo_head_data(fname, halofinder="rockstar", feature="mvir", zspace=False)->Table:
    if halofinder.lower() == "rockstar":
        return load_rockstar_head_data(fname, feature, zspace)
    
def load_rockstar_head_data(fname:str, feature:str, zspace:bool=False, los:str="z")->Table:
    fname = os.path.join(fname, "out_0.list")
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"File {fname} does not exist!")

    f = open(fname, "r")
    header = []
    while(1):
        line = f.readline()
        if line[0] == "#":
            header.append(line[:-1])
        else:
            break
    boxsize = float(header[6].split(" ")[-2])
    scale_factor = float(header[1].split("=")[1])
    redshift = 1./scale_factor - 1
    OmegaM = float(header[2].split(";")[0].split("=")[1])

    halo = Table(meta={'boxsize': boxsize, 'redshift': redshift, 'OmegaM': OmegaM})
    tmp = np.loadtxt(fname)

    halo["x"] = tmp[:,rockstar_idx_dict["x"]]
    halo["y"] = tmp[:,rockstar_idx_dict["y"]]
    halo["z"] = tmp[:,rockstar_idx_dict["z"]]

    if zspace:
        s_los = cal_zrsd(OmegaM, redshift, halo[:,los], tmp[:,rockstar_idx_dict[los]])
        halo[los] = s_los

    halo[feature] = tmp[:,rockstar_idx_dict[feature]]
    
    return halo

def SHAM_sigma_model(sigma, pos, ftr_val, boxsize, ref_num_den=3.5e-4, seed=None):
    rng = np.random.default_rng(seed=seed)

    Nhalo = len(pos)
    scatter = rng.normal(loc=0.0, scale=sigma, size=Nhalo)
    idx_plus = (scatter > 0)
    idx_minus = (scatter < 0)
    scatter[idx_plus] = scatter[idx_plus] + 1
    scatter[idx_minus] = np.exp(scatter[idx_minus])

    Ntarget = int(boxsize*boxsize*boxsize*ref_num_den)

    ftr_scat = ftr_val*scatter

    idxed_arr = np.c_[pos, ftr_scat]
    gsamples = heapq.nlargest(Ntarget, idxed_arr, key=lambda x:x[-1])
    gsamples = np.asarray(gsamples)

    return gsamples