import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec

rdir = "results/Vprof/cosmo1/a_0.7692/SHAM0/"
Rmins = np.linspace(15, 24, 10)
Rmaxs = np.linspace(16, 25, 10)
# for ibin in range(10):
#     Rvi = 0.5*(Rmins[ibin]+Rmaxs[ibin])
#     prof_std = np.loadtxt(os.path.join(rdir, f"rvbin{ibin}_std.txt"))[:,1]
#     prof_1per = np.loadtxt(os.path.join(rdir, f"rvbin{ibin}_1percent.txt"))[:,1]
#     prof_mesh = np.loadtxt(os.path.join(rdir, f"rvbin{ibin}_mesh_cat.txt"))[:,1]
#     Rv = np.logspace(np.log10(0.1*Rvi), np.log10(3*Rvi), len(prof_std))

#     fig = plt.figure()
#     gs = gridspec.GridSpec(5,5)
#     ax1 = fig.add_subplot(gs[:-1,:])
#     ax1.plot(Rv, prof_std, "--k", label="std")
#     ax1.plot(Rv, prof_1per, "r", label=r"1% sample")
#     ax1.plot(Rv, prof_mesh, "b", label="Nmesh=512")
#     ax1.legend()
#     ax1.set_xscale("log")
#     ax1.set_xticklabels([])
#     ax1.set_ylabel(r"$\delta$")
#     ax2 = fig.add_subplot(gs[-1,:])
#     ax2.plot(Rv, (1+prof_1per)/(1+prof_std)-1, "r")
#     ax2.plot(Rv, (1+prof_mesh)/(1+prof_std)-1, "b")
#     ax2.plot([Rv[0], Rv[-1]], [0,0], "--k")
#     ax2.set_xscale("log")
#     ax2.set_xlabel("Rv")
#     ax2.set_ylabel(r"$\Delta$")
#     plt.savefig(f"compare_rvbin{ibin}.png", dpi=400)

prof_std = np.loadtxt(os.path.join(rdir, f"stacked_std.txt"))[:,1]
prof_1per = np.loadtxt(os.path.join(rdir, f"stacked_1percent.txt"))[:,1]
prof_mesh = np.loadtxt(os.path.join(rdir, f"stacked_mesh_cat.txt"))[:,1]
Rv = np.logspace(-1, np.log10(3), len(prof_std))

fig = plt.figure()
gs = gridspec.GridSpec(5,5)
ax1 = fig.add_subplot(gs[:-1,:])
ax1.plot(Rv, prof_std, "--k", label="std")
ax1.plot(Rv, prof_1per, "r", label=r"1% sample")
# ax1.plot(Rv, prof_mesh, "b", label="Nmesh=512")
ax1.legend()
ax1.set_xscale("log")
ax1.set_xticklabels([])
ax1.set_ylabel(r"$\delta$")
ax2 = fig.add_subplot(gs[-1,:])
ax2.plot(Rv, (1+prof_1per)/(1+prof_std)-1, "r")
# ax2.plot(Rv, (1+prof_mesh)/(1+prof_std)-1, "b")
ax2.plot([Rv[0], Rv[-1]], [0,0], "--k")
ax2.set_xlabel("Rv")
ax2.set_xscale("log")
ax2.set_ylabel(r"$\Delta$")
plt.savefig(f"compare_stacked.png", dpi=400)