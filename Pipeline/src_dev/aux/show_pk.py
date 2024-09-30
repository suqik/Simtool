import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
# from nbodykit.lab import FFTPower
from scipy.interpolate import interp1d

# fdir = "/public/share/ace66so15x/suchen/L1000_N1024_rlzs/"

# for i in range(10):
#     tmp = np.loadtxt(fdir+f"rlz{i}/powerspec/powerspec_1.0000.txt")
#     plt.loglog(tmp[1:,0], tmp[1:,1])

# plt.savefig("power_rlzs.png", dpi=400)

pk_fid = np.loadtxt("Pk_measurement_fiducial_fid.txt")
pk_inv = np.loadtxt("Pk_measurement_fiducial_inv.txt")
pk_th  = np.loadtxt("Pk_fiducial_theory.txt")
pk_interp = interp1d(pk_th[:,0], pk_th[:,1])

fig = plt.figure()
gs = GridSpec(5,5)
ax1 = fig.add_subplot(gs[:-1,:])
ax1.loglog(pk_fid[:,0], pk_fid[:,1], lw=0.5)
ax1.loglog(pk_inv[:,0], pk_inv[:,1], lw=0.5)
ax1.loglog(pk_fid[:,0], pk_interp(pk_fid[:,0]), lw=0.5, c='k')
ax1.set_xticklabels([])
ax1.set_ylim(80, 3e4)

ax2 = fig.add_subplot(gs[-1,:])
ax2.semilogx(pk_fid[:,0], pk_fid[:,1]/pk_interp(pk_fid[:,0])-1)
ax2.semilogx(pk_inv[:,0], pk_inv[:,1]/pk_interp(pk_fid[:,0])-1)
ax2.semilogx([pk_fid[0,0], pk_fid[-1,0]], [0,0], '--k')
ax2.set_ylim(-0.1,0.1)

fig.savefig("power_fix_pair_test_fix.png", dpi=400)

# results_fid = FFTPower.load("Pk_IC_fiducial_fid.json")
# results_inv = FFTPower.load("Pk_IC_fiducial_inv.json")

# plt.loglog(results_fid.power["k"], results_fid.power["power"])
# plt.loglog(results_inv.power["k"], results_inv.power["power"], "--")
# plt.savefig("/public/home/suchen/Programs/Simtool/Pipeline/figs/gen_IC/test_pk.png", dpi=400)