import sys
sys.path.append("/public/home/suchen/Programs/Simtool/Pipeline/src_dev/")
import numpy as np
from nbodykit.source.catalog import BigFileCatalog
from nbodykit.cosmology import HalofitPower, Cosmology
from utils.power import get_powerspec_cat

outdir = "/public/home/suchen/Programs/Simtool/Pipeline/src_dev/aux/"

# fcosmo_list = "/public/home/suchen/Programs/Simtool/Pipeline/cfgs/validation/val_cosmo_list.txt"
# cosmo_list = {}
# with open(fcosmo_list, "r") as f:
#     line = f.readline().strip("\n")
#     key_val = line.split(" ")[1:]
#     for ikv in key_val:
#         key, val =  ikv.split("=")
#         cosmo_list[key] = float(val)
#     line = f.readline().strip("\n")
#     vary_param_name = line.split(" ")[1:]

# vary_cosmo_val = np.loadtxt(fcosmo_list)
# if len(vary_param_name) == 1:
#     cosmo_list[vary_param_name] = vary_cosmo_val
# else:
#     for idx, iname in enumerate(vary_param_name):
#         cosmo_list[iname] = vary_cosmo_val[:,idx]

# validate_label = 50
### theory
def mk_fin_Pk(hubble, Om0, Ob0, ns, sigma8, output):
    MYcosmo = Cosmology(m_ncdm=[],
                        Omega0_b=Ob0,
                        Omega0_cdm=Om0 - Ob0,
                        h=hubble)\
                        .match(sigma8=sigma8)

    pklin = HalofitPower(MYcosmo, redshift=0.3)
    k = np.logspace(-3,1,1000, endpoint=True)
    np.savetxt(output, np.c_[k, pklin(k)])

    return None

# mk_fin_Pk(cosmo_list["hubble"], 
#           cosmo_list["OmegaM"],
#           cosmo_list["Omegab"], 
#           cosmo_list["ns"], 
#           cosmo_list["S8"]/np.sqrt(cosmo_list["OmegaM"]/0.3),
#           f"Pk_theory_cosmo{validate_label}.txt")

# mk_fin_Pk(0.6727, 0.3156, 0.0491, 0.9667, 0.831, 
#           outdir+"Pk_fiducial_theory.txt")

### measurement
# catdir = "/public/share/ace66so15x/suchen/L1000_N1024_validation/"
# data = BigFileCatalog(catdir+f"cosmo{validate_label}/a_0.7692/", dataset="1/", header="Header")
# kout, Pkout = get_powerspec_cat(data["Position"].compute(), Nmesh=512, 
#                                 boxsize=1000., threads=32)

# np.savetxt(f"Pk_measurement_cosmo{validate_label}.txt", np.c_[kout, Pkout[0]])

catdir = "/public/share/ace66so15x/suchen/L1000_N1024_fix_pair/"
data = BigFileCatalog(catdir+f"fiducial/fid/a_0.7692/", dataset="1/", header="Header")
kout, Pkout = get_powerspec_cat(data["Position"].compute(), Nmesh=512, 
                                boxsize=1000., threads=32)

np.savetxt(outdir+"Pk_measurement_fiducial_fid.txt", np.c_[kout, Pkout[0]])

catdir = "/public/share/ace66so15x/suchen/L1000_N1024_fix_pair/"
data = BigFileCatalog(catdir+f"fiducial/inv/a_0.7692/", dataset="1/", header="Header")
kout, Pkout = get_powerspec_cat(data["Position"].compute(), Nmesh=512, 
                                boxsize=1000., threads=32)

np.savetxt(outdir+"Pk_measurement_fiducial_inv.txt", np.c_[kout, Pkout[0]])