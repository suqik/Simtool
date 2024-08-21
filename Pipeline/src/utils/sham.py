import heapq
import numpy as np
from nbodykit.source.catalog import BigFileCatalog

class halo_info(object):
    def __init__(self, header, x, y, z, vx=None, vy=None, vz=None, mass=None, Vmax=None, pid=None):
        self.header = header
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.mass = mass
        self.Vmax = Vmax
        self.pid = pid
        self.size = len(self.x)

    def confirm_position(self, boxsize=None):
        if 'boxsize' in self.header.keys():
            if boxsize is not None and np.abs(boxsize - self.header['boxsize']) < 1e-5:
                print("Input boxsize does not match the header! Will ignore input value.")
            boxsize = self.header['boxsize']
        else:
            if boxsize is None:
                raise ValueError("Need boxsize as input!")
        
        self.x = (self.x + boxsize)%boxsize
        self.y = (self.y + boxsize)%boxsize
        self.z = (self.z + boxsize)%boxsize
        if hasattr(self, 'xrsd'):
            self.xrsd = (self.xrsd + boxsize)%boxsize
        if hasattr(self, 'yrsd'):
            self.yrsd = (self.yrsd + boxsize)%boxsize
        if hasattr(self, 'zrsd'):
            self.zrsd = (self.zrsd + boxsize)%boxsize
    
    def cal_zrsd(self, OmegaM=None, redshift=None, los='z'):
        if 'OmegaM' in self.header.keys():
            if OmegaM is not None and OmegaM != self.header['OmegaM']:
                print("Input OmegaM does not match the header! Will ignore input value.")
            OmegaM = self.header['OmegaM']
        else:
            if OmegaM is None:
                raise ValueError("Need OmegaM as input!")

        if 'redshift' in self.header.keys():
            if redshift is not None and redshift != self.header['redshift']:
                print("Input redshift does not match the header! Will ignore input value.")
            redshift = self.header['redshift']
        else:
            if redshift is None:
                raise ValueError("Need redshift as input!")
        
        hubble = 100 * np.sqrt(OmegaM * (redshift + 1.0)**3 + (1-OmegaM))
        rsd_shift = (redshift + 1.0) / hubble

        if los == 'x':
            self.xrsd = self.x + self.vx*rsd_shift
        if los == 'y':
            self.yrsd = self.y + self.vy*rsd_shift
        if los == 'z':
            self.zrsd = self.z + self.vz*rsd_shift

def load_bigfile_halo(idir):
    tmp = BigFileCatalog(idir, dataset="LL-0.200/", header="Header")
    attrs = tmp.attrs
    particle_mass = tmp.attrs['M0'][0]*1e10

    mass = (tmp['Length']*particle_mass)
    header = {
        'boxsize': attrs['BoxSize'][0],
        'redshift': 1./(attrs['ScalingFactor'])[0] - 1,
        'OmegaM': attrs['OmegaM'][0],
        'PMass': attrs['M0'][0]*1e10
    }
    halo = halo_info(
        header=header,
        x=tmp['Position'][:,0].compute(),
        y=tmp['Position'][:,1].compute(),
        z=tmp['Position'][:,2].compute(),
        vx=tmp['Velocity'][:,0].compute(),
        vy=tmp['Velocity'][:,1].compute(),
        vz=tmp['Velocity'][:,2].compute(),
        mass=mass.compute(),
        pid=np.ones(tmp.csize)-1
    )
    halo.cal_zrsd()
    halo.confirm_position()

    return halo

def load_rockstar_halo(fname, feature="Mass", z_space=True):
    f = open(fname, "r")
    header = []
    while(1):
        line = f.readline()
        if line[0] == "#":
            header.append(line[:-1])
        else:
            break
    
    scale_factor = float(header[1].split("=")[1])
    redshift = 1./scale_factor - 1
    OmegaM = float(header[2].split(";")[0].split("=")[1])
    boxsize = float(header[6].split(" ")[-2])

    header = {
        'boxsize': boxsize,
        'redshift': redshift,
        'OmegaM': OmegaM,
    }
    
    if feature == "Mass":
        halo_mass = []
    if feature == "Vmas":
        halo_Vmax = []
    halo_x = []
    halo_y = []
    halo_z = []
    halo_vz = []
    if feature == "Mass":
        for line in f.readlines():
            line = line.split(" ")
            halo_mass.append(float(line[2]))
            halo_x.append(float(line[8]))
            halo_y.append(float(line[9]))
            halo_z.append(float(line[10]))
            halo_vz.append(float(line[13]))
        
        halo_mass = np.asarray(halo_mass)
        halo_x = np.asarray(halo_x)
        halo_y = np.asarray(halo_y)
        halo_z = np.asarray(halo_z)
        halo_vz = np.asarray(halo_vz)
        halo = halo_info(
            header=header,
            x=halo_x,
            y=halo_y,
            z=halo_z,
            vz=halo_vz,
            mass=halo_mass
        )

    if feature == "Vmax":
        for line in f.readlines():
            line = line.split(" ")
            halo_Vmax.append(float(line[3]))
            halo_x.append(float(line[8]))
            halo_y.append(float(line[9]))
            halo_z.append(float(line[10]))
            halo_vz.append(float(line[13]))
        
        halo_Vmax = np.asarray(halo_Vmax)
        halo_x = np.asarray(halo_x)
        halo_y = np.asarray(halo_y)
        halo_z = np.asarray(halo_z)
        halo_vz = np.asarray(halo_vz)

        halo = halo_info(
            header=header,
            x=halo_x,
            y=halo_y,
            z=halo_z,
            vz=halo_vz,
            Vmax=halo_Vmax
        )
    
    if z_space:
        halo.cal_zrsd()
    halo.confirm_position()

    return halo

def SHAM_model(sigma, catalog:halo_info, feature="Mass", z_space=True, ref_num_den=3.5e-4, seed=None, rng=None):
    if rng is None:
        rng = np.random.default_rng(seed=seed)

    scatter = rng.normal(loc=0.0, scale=sigma, size=catalog.size)
    idx_plus = (scatter > 0)
    idx_minus = (scatter < 0)
    scatter[idx_plus] = scatter[idx_plus] + 1
    scatter[idx_minus] = np.exp(scatter[idx_minus])

    boxsize = catalog.header['boxsize']
    Ntarget = int(boxsize*boxsize*boxsize*ref_num_den)

    if feature == "Mass":
        ftr = catalog.mass*scatter
    if feature == "Vmax":
        ftr = catalog.Vmax*scatter

    if z_space:
        idxed_arr = np.c_[catalog.x, catalog.y, catalog.zrsd, ftr]
    else:
        idxed_arr = np.c_[catalog.x, catalog.y, catalog.z, ftr]
    gsamples = heapq.nlargest(Ntarget, idxed_arr, key=lambda x:x[-1])
    gsamples = np.asarray(gsamples)

    return gsamples