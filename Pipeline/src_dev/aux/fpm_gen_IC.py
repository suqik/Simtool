import sys
sys.path.append("/public/home/suchen/Programs/Simtool/Pipeline/src_dev/")
from utils.mk_conf_func import mk_ini_lineark
import datetime

start = datetime.datetime.now()

### generate fix and pair initial density field
seed = 7296
boxsize = 1000
Nmesh = 1024
fid_mesh_r = mk_ini_lineark(hubble=0.6727, 
               Om0=0.3156, 
               Ob0=0.0491, 
               ns=0.9667, 
               sigma8=0.831, 
               Nmesh=Nmesh, 
               boxsize=boxsize,
               output="/public/share/ace66so15x/suchen/L1000_N1024_fix_pair/fiducial/IC/fid", 
               unitary_amplitude=True, 
               seed=seed,)

inv_mesh_r = mk_ini_lineark(hubble=0.6727, 
               Om0=0.3156, 
               Ob0=0.0491, 
               ns=0.9667, 
               sigma8=0.831, 
               Nmesh=Nmesh, 
               boxsize=boxsize,
               output="/public/share/ace66so15x/suchen/L1000_N1024_fix_pair/fiducial/IC/inv",
               unitary_amplitude=True, 
               inverted_phase=True, 
               seed=seed,)

end = datetime.datetime.now()
print(end - start)

### test the saved field
# seed = 7296
# boxsize = 100
# Nmesh = 64
# fid_mesh_r = mk_ini_lineark(hubble=0.6727, 
#                Om0=0.3156, 
#                Ob0=0.0491, 
#                ns=0.9667, 
#                sigma8=0.831, 
#                Nmesh=Nmesh, 
#                boxsize=boxsize,
#                output="fid/", 
#                unitary_amplitude=True,
#                seed=seed,)

# inv_mesh_r = mk_ini_lineark(hubble=0.6727, 
#                Om0=0.3156, 
#                Ob0=0.0491, 
#                ns=0.9667, 
#                sigma8=0.831, 
#                Nmesh=Nmesh, 
#                boxsize=boxsize,
#                output="inv/",
#                unitary_amplitude=True,
#                inverted_phase=True, 
#                seed=1234,)

# end = datetime.datetime.now()
# print(end - start)

# from matplotlib import pyplot as plt
# from nbodykit.lab import BigFileMesh, FFTPower

# results = FFTPower(fid_mesh_r, mode="1d")
# results.save("Pk_IC_fiducial_fid.json")
# results = FFTPower(inv_mesh_r, mode="1d")
# results.save("Pk_IC_fiducial_inv.json")

# fig = plt.figure()
# ax1 = fig.add_subplot(121)
# ax1.imshow(fid_mesh_r.preview(Nmesh=256, axes=(0,1)))
# ax2 = fig.add_subplot(122)
# ax2.imshow(inv_mesh_r.preview(Nmesh=256, axes=(0,1)))

# fig.savefig("/public/home/suchen/Programs/Simtool/Pipeline/figs/gen_IC/test.png", dpi=400)
