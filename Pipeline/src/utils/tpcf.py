import numpy as np
from pyfcfc.boxes import py_compute_cf
from pyfcfc.utils import add_pair_counts

def get_auto_wp(data, boxsize, bins_info, pyfcfc_conf, rand=None, seed=None, wdat=None, wran=None):
    dtype = data.dtype
    dsize = len(data)
    if wdat is None:
        wdat = np.ones(dsize).astype(dtype)
    if rand is None:
        rng = np.random.default_rng(seed=seed)
        rand = rng.random((dsize*10, 3), dtype=dtype)*boxsize
        wran = np.ones(dsize*10).astype(dtype)

    min_sep = bins_info["min_sep"]
    max_sep = bins_info["max_sep"]
    nbins   = bins_info["nbins"]
    if "min_pi" in bins_info.keys():
        min_pi  = bins_info["min_pi"]
    else:
        min_pi  = 0
    if "max_pi" in bins_info.keys():
        max_pi  = bins_info["max_pi"]
    else:
        max_pi  = 40
    if "npbins" in bins_info.keys():
        npbins  = bins_info["npbins"]
    else:
        npbins  = 40

    results = py_compute_cf([data, rand], [wdat, wran],
                            np.logspace(np.log10(min_sep), np.log10(max_sep), int(nbins)), 
                            np.linspace(min_pi, max_pi, int(npbins)), 
                            1,
                            label=["D", "R"],
                            box=boxsize,
                            pair=["DD", "DR", "RR"],
                            cf = "(DD - 2*DR + RR)/RR",
                            conf = pyfcfc_conf
                            )

    return results

def get_cross_iso(cat1, cat2, boxsize,
                  wcat1=None, wcat2=None,
                    bins_info = {
                    "min_sep": 5,
                    "max_sep": 35,
                    "nbins": 15,
                    }
                  ):
    
    nbins = bins_info['nbins']
    min_sep = bins_info['min_sep']
    max_sep = bins_info['max_sep']

    dtype = cat1.dtype
    if wcat1 is None:
        wcat1 = np.ones(len(cat1)).astype(dtype)
    else:
        wcat1 = wcat1.astype(dtype)

    cat2 = cat2.astype(dtype)
    if wcat2 is None:
        wcat2 = np.ones(len(cat2)).astype(dtype)
    else:
        ### FIXME Check this !!! (dimension)
        wcat2 = wcat2.astype(dtype)

    ### measure iso 2pcf
    results = py_compute_cf([cat1, cat2], [wcat1, wcat2], 
                            10**(np.linspace(np.log10(min_sep), np.log10(max_sep), nbins)), 
                            None, 
                            100, 
                            label = ['D', 'V'], # Catalog labels matching the number of catalogs provided
                            bin=1, # bin type for multipoles
                            pair = ['DV'], # Desired pair counts
                            box=boxsize, 
                            multipole = [0], # Multipoles to compute
                            cf = ['DV / @@ - 1'], # CF estimator (not necessary if only pair counts are required)
                            verbose = 'F'
                            ) 
    
    return results

def get_stacked_vprof(dm, void, boxsize, Rmins, Rmaxs, Rv_col=-1, wdm=None, wv=None,
                      min_sep_inRv = 0.1, max_sep_inRv=3, nbins=15,
                      ind_file=None, stack_file=None):
    
    Ndm = len(dm)

    ### choose void catalog in considered size bins
    general_cut = (void[:,Rv_col] > Rmins[0]) & (void[:,Rv_col] < Rmaxs[-1])
    void = void[general_cut]
    del general_cut

    results_list = {"s_list":[],
                    "xi_iso_list":[]}
    
    for iR in range(len(Rmins)):
        Rl = Rmins[iR]
        Ru = Rmaxs[iR]
        Rv = (Ru + Rl)/2.

        min_sep = min_sep_inRv*Rv
        max_sep = max_sep_inRv*Rv

        dlgbin = np.log10(max_sep/min_sep)/(nbins-1)
        min_sep = 10**(np.log10(min_sep)-0.5*dlgbin)
        max_sep = 10**(np.log10(max_sep)+0.5*dlgbin)
        nedges = nbins + 1

        bins_info = {
        "min_sep": min_sep,
        "max_sep": max_sep,
        "nbins": nedges,
        }
        
        vcut = (void[:,Rv_col]>Rl) & (void[:,Rv_col]<Ru)
        Nv = np.sum(vcut)

        results = get_cross_iso(dm, void[vcut][:,:Rv_col], boxsize, 
                                wcat1=wdm, wcat2=wv, 
                                bins_info=bins_info)
        
        sep = results['s']
        sep_inRv = sep/Rv
        xi_iso = np.array(results['multipoles'])[0,0,:]
        results_list["s_list"].append(sep)
        results_list["xi_iso_list"].append(xi_iso)
        # if ind_file is not None:
        #     f = open(ind_file+f"rvbin{iR}.txt", "w+", encoding="utf-8")
        #     f.write("# scen (Mpc/h) scen (Rv) xi_iso (Rvmin={:.2f} Rvmax={:.2f} (Mpc/h))\n".format(Rl, Ru))
        #     np.savetxt(f, np.c_[sep, sep/Rv, xi_iso])
        #     f.close()
        
        edges = np.logspace(np.log10(min_sep), np.log10(max_sep), nedges)
        RR_ana_k = 4./3.*np.pi*((edges[1:]/boxsize)**3-(edges[:-1]/boxsize)**3)*Nv*Ndm
        RR_ana = RR_ana + RR_ana_k if iR > 0 else RR_ana_k
        
        total_results = add_pair_counts(total_results, results) if iR > 0 else results

    xi_iso_stacked = np.sum(total_results['pairs']['DV'], axis=1)*total_results['normalization']['DV']/RR_ana - 1
    results_list["sep_inRv"] = sep_inRv
    results_list["xi_iso_stacked"] = xi_iso_stacked
    # if stack_file is not None:
    #     f = open(stack_file+"stacked.txt", "w+", encoding="utf-8")
    #     f.write("# scen (Rv) xi_iso (Rvmin={:.2f} Rvmax={:.2f} (Mpc/h))\n".format(Rmins[0], Rmaxs[-1]))
    #     np.savetxt(f, np.c_[sep_inRv, xi_iso])
    #     f.close()
    
    return results_list
