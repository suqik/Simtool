import os, sys
import copy
import configparser
import numpy as np
sys.path.append("/public/home/suchen/Programs/Simtools/Pipeline/src_dev/")
import random
from nbodykit.source.catalog import BigFileCatalog

from utils.mk_conf_func import *
from utils.convert import *
from utils.sham import *
from utils.tpcf import *
from utils.power import *

'''
TODO: Reconstruct this part
1. Each runner only processes ONE cosmology, all of the cosmology routines 
should be move to external scripts (like run_fastpm.py, run_gal_void.py, etc.).
2. For FastPM, just reading cosmological parameters from FILES, consider if 
support accept parameters direct from functional parameter.
3. For cosmological realizations and SHAM parameter/realizations, consider them
carefully later. (Maybe realizations as internal functions, parameters external?)
To Be Continue -->
'''

class Base_Runner(object):
    def __init__(self, conf=None):
        if conf is None:
            self.conf = configparser.ConfigParser()
        else:
            self.conf = conf
    
    def load_config_file(self, conf_file):
        if not os.path.isfile(os.path.join(os.getcwd(), conf_file)):
            print("Configure file does not exist!")
            exit()
        
        self.conf = configparser.ConfigParser()
        self.conf.read(conf_file)

class PREPARE_Runner(Base_Runner):
    def __init__(self, conf=None):
        super(PREPARE_Runner, self).__init__(conf)
    
    def set_params(self, **kwargs):
        ### general params
        self.cfgbase = self.conf.get("General", "cfgbase").strip("\"")
        self.cfgsubbase = self.conf.get("General", "cfgsubbase").strip("\"")
        
        ### cosmology-related params
        if "input" in self.conf["General"]:
            self.FROMFILE = True
            self.input = self.conf.get("General", "input").strip("\"")
        else:
            self.FROMFILE = False
            self.nparams = self.conf["General"].getint("nparams")
            self.ncosmo = self.conf["General"].getint("ncosmo")
            self.fix_cosmo_names = list(map(str, self.conf.get("General", "fix_cosmo_names").split(", ")))
            self.fix_cosmo_vals  = list(map(float, self.conf.get("General", "fix_cosmo_vals").split(", ")))            
            self.param_names = list(map(str, self.conf.get("General", "param_names").split(", ")))
            self.prior_low = list(map(float, self.conf.get("General", "prior_low").split(", ")))
            self.prior_up  = list(map(float, self.conf.get("General", "prior_up").split(", ")))
            self.cosmo_param_seed = self.conf["General"].getint("seed")
            self.output = self.conf.get("General", "output").strip("\"")

        if "fix_seed" in self.conf["General"]:
            self.fix_seed = self.conf["General"].getboolean("fix_seed")
        else:
            self.fix_seed = False

        ### halo-related params
        self.ROCKSTAR = self.conf["FastPM"].getboolean("ROCKSTAR")
        self.cvt_nfile  = self.conf["Convert"].getint("nfile")
        self.cvt_opbase = self.conf.get("Convert", "outputbase").strip("\"")

        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

    def get_cosmo_params(self):
        if self.FROMFILE:
            cosmo_file_path = self.cfgbase+self.input
            if os.path.isfile(cosmo_file_path):
                self.fix_cosmo_names = []
                self.fix_cosmo_vals = []
                with open(cosmo_file_path, "r") as f:
                    tmp = f.readline().strip("\n").split(" ")[1:]
                    for key_val in tmp:
                        key = key_val.split("=")[0]
                        val = float(key_val.split("=")[1])
                        self.fix_cosmo_names.append(key)
                        self.fix_cosmo_vals.append(val)
                    tmp = f.readline().strip("\n").split(" ")[1:]
                    self.param_names = tmp
                cosmo_params = np.loadtxt(cosmo_file_path)
            if len(cosmo_params.shape) == 1:
                cosmo_params = cosmo_params[np.newaxis,:]
            ip_ncosmo, ip_nparams = cosmo_params.shape
            self.ncosmo = ip_ncosmo
            self.nparams = ip_nparams
        else:
            ### execute
            cosmo_param_rng = np.random.default_rng(seed=self.cosmo_param_seed)
            cosmo_params = cosmo_param_rng.uniform(
                low=self.prior_low, 
                high=self.prior_up, 
                size=(self.ncosmo, self.nparams)
                )

            f = open(self.cfgbase+self.output, "w+", encoding="utf-8")
            f.write("# {}\n".format(" ".join([f"{self.fix_cosmo_names[i]}={self.fix_cosmo_vals[i]:.4f}" for i in range(len(self.fix_cosmo_names))])))
            f.write("# {}\n".format(" ".join(self.param_names)))
            np.savetxt(f, cosmo_params, fmt="%3f %3f")
            f.close()

        param_dict = {}
        for i in range(self.nparams):
            param_dict[self.param_names[i]] = cosmo_params[:,i]

        return copy.deepcopy(param_dict)
    
    def declare(self):
        ### TODO: write declare
        print("Beginning preparing necessary configuration files.", flush=True)
        if self.FROMFILE:
            print(f"Load cosmological parameters from file {self.cosmo_file}", flush=True)
            print("Fixed parameters: {}".format(" ".join([f"{self.fix_cosmo_names[i]}={self.fix_cosmo_vals[i]:.4f}" for i in range(len(self.fix_cosmo_names))])), flush=True)
            print(f"Varied parameters: {self.param_names}", flush=True)
        else:
            print(f"Randomly sample parameters and will save to file {self.output}", flush=True)
            print("Fixed parameters: {}".format(" ".join([f"{self.fix_cosmo_names[i]}={self.fix_cosmo_vals[i]:.4f}" for i in range(len(self.fix_cosmo_names))])), flush=True)
            print(f"Varied parameters: {self.param_names}", flush=True)

        print(f"fix_seed: {self.fix_seed}")

    def run(self, relic=""):
        vary_dict = self.get_cosmo_params()
        relic_input = str(relic)

        param_dict = {}
        for i in range(len(self.fix_cosmo_names)):
            param_dict[self.fix_cosmo_names[i]] = self.fix_cosmo_vals[i]

        for icosmo in range(self.ncosmo):
            if relic_input == "" and self.ncosmo > 1:
                relic = f"{icosmo}"
            # print(relic)
            ### Generate FastPM configuration file ###
            fpm_cfgpath = self.cfgbase+self.cfgsubbase+f"{relic}/fastpm/"
            if not os.path.isdir(fpm_cfgpath):
                os.makedirs(fpm_cfgpath)

            # update cosmology parameter
            for iname in range(self.nparams):
                param_dict[self.param_names[iname]] = vary_dict[self.param_names[iname]][icosmo]

            hubble = param_dict["hubble"]
            Om0 = param_dict["OmegaM"]
            Ob0 = param_dict["Omegab"]
            ns  = param_dict["ns"]
            S8  = param_dict["S8"]
            # calculate initial matter power spectrum        
            sigma8 = S8/np.sqrt(Om0/0.3)

            mk_ini_Pk(hubble, Om0, Ob0, ns, sigma8, fpm_cfgpath+"Pkini.txt")

            # write fastpm conf    
            if self.fix_seed:        
                mk_fastpm_conf(self.conf, 0, Om0, hubble, fpm_cfgpath+"Pkini.txt", fpm_cfgpath+"fpm.lua", relic=relic)
            else:
                mk_fastpm_conf(self.conf, icosmo, Om0, hubble, fpm_cfgpath+"Pkini.txt", fpm_cfgpath+"fpm.lua", relic=relic)

            if self.ROCKSTAR:
                rstar_cfgpath = self.cfgbase+self.cfgsubbase+f"{relic}/rockstar/"
                if not os.path.isdir(rstar_cfgpath):
                    os.makedirs(rstar_cfgpath)
                redshifts = list(map(float, self.conf.get("FastPM", "redshifts").split(", ")))
                for idx, zi in enumerate(redshifts):
                    mk_rockstar_conf(self.conf, 
                                    cvt_opbase=self.cvt_opbase.split("/")[0], 
                                    cvt_nfile=self.cvt_nfile, 
                                    filename=self.cvt_opbase.split("/")[1]+"{:03d}.<block>".format(idx), 
                                    output=rstar_cfgpath+"z{:.2f}.cfg".format(zi), 
                                    relic=relic, 
                                    redshift=zi)
            
class FASTPM_Runner(Base_Runner):
    def __init__(self, conf=None):
        super(FASTPM_Runner, self).__init__(conf)

    def set_params(self, **kwargs):
        self.base = self.conf.get("General","cfgbase").strip("\"")
        self.subbase = self.conf.get("General","cfgsubbase").strip("\"")
        if "input" in self.conf["General"]:
            fcosmo_param = self.conf.get("General", "input").strip("\"")
            tmp = np.loadtxt(self.base+fcosmo_param)
            if len(tmp.shape) == 1:
                self.ncosmo = 1
            else:
                self.ncosmo = tmp.shape[0]
        else:
            self.ncosmo = self.conf["General"].getint("ncosmo")

        self.FastPM_exec = self.conf.get("FastPM", "FastPM_EXEC").strip("\"")

        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

    def run(self, nCPUs=None, snapname_relic="", iteration=False):
        if nCPUs is None or nCPUs == 1:
            cmdbase = self.FastPM_exec+" "
        else:
            cmdbase = "mpirun -np "+f"{nCPUs} "+self.FastPM_exec+" "
        os.environ["OMP_NUM_THREADS"] = "1"
        if not iteration:
            fpm_cfgpath = self.base+self.subbase+f"{snapname_relic}/fastpm/fpm.lua"
            if not os.path.isfile(fpm_cfgpath):
                exit()
            cmd = cmdbase + fpm_cfgpath
            print(cmd, flush=True)
            os.system(cmd)
        else:
            for icosmo in range(self.ncosmo):
                snapname_relic = f"{icosmo}"
                fpm_cfgpath = self.base+self.subbase+f"{snapname_relic}/fastpm/fpm.lua"
                cmd = cmdbase + fpm_cfgpath
                print(cmd, flush=True)
                os.system(cmd)

class ROCKSTAR_Runner(Base_Runner):
    def __init__(self, conf=None):
        super(ROCKSTAR_Runner, self).__init__(conf)

    def set_params(self, **kwargs):
        ### Convert params
        self.nfile = self.conf["Convert"].getint("nfile")
        self.precision = self.conf.get("Convert", "precision").strip("\"")
        self.gadgetbase = self.conf.get("Convert", "outputbase").strip("\"")

        ### Input & output
        self.ipdir = self.conf.get("FastPM", "snapdir").strip("\"")
        self.ipbase = self.conf.get("FastPM", "snapbase").strip("\"")
        self.redshifts = list(map(float, self.conf.get("FastPM", "redshifts").split(", ")))

        ### Rockstar part
        self.boxsize = self.conf["FastPM"].getfloat("boxsize")
        self.base = self.conf.get("General","cfgbase").strip("\"")
        self.subbase = self.conf.get("General","cfgsubbase").strip("\"")
        self.rstaropbase = self.conf.get("ROCKSTAR", "outputbase").strip("\"")

        self.Rockstar_exec = self.conf.get("ROCKSTAR", "ROCKSTAR_EXEC").strip("\"")

        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

    def run(self, snapname_relic=""):
        for irds, zi in enumerate(self.redshifts):
            snappath = self.ipdir+self.ipbase+"{}/a_{:.4f}/".format(snapname_relic, 1./(1+zi))
            gadgetpath = snappath+self.gadgetbase+"{:03d}".format(irds)

            ### execute convert
            if not os.path.isfile(gadgetpath+".0"):
                Convert(snappath, gadgetpath, self.nfile, self.precision)

            ### run rockstar
            rstar_cfgpath = self.base+self.subbase+"{}/rockstar/z{:.2f}.cfg".format(snapname_relic, zi)
            if not os.path.isdir(snappath+self.rstaropbase):
                os.mkdir(snappath+self.rstaropbase)

            ### I don't know why but this way works ... 
            tmp_fname = "tmp_run_rstar_{}.sh".format(snapname_relic)
            f = open(tmp_fname, "w+")
            
            f.write("#!/bin/bash\n")
            f.write(f"RSTAR={self.Rockstar_exec}\n")
            # f.write(f"FINDPAR={FindPAR_exec}\n")
            f.write(f"BOXSIZE={self.boxsize}\n")
            f.write(f"CFG={rstar_cfgpath}\n")
            f.write(f"IDIR={gadgetpath}\n")
            f.write(f"ODIR={snappath+self.rstaropbase}\n")
            f.write("$RSTAR -c $CFG &\n")
            f.write("export ODIR\n")
            f.write("perl -e \'sleep 1 while (!(-e \"$ENV{ODIR}/auto-rockstar.cfg\"))\'\n")
            f.write("$RSTAR -c ${ODIR}/auto-rockstar.cfg\n")
            # f.write("$FINDPAR ${ODIR}/out_0.list $BOXSIZE >${ODIR}/out_0_wsub.list\n")
            f.write("rm -r ${ODIR}/halos* ${ODIR}/*.cfg ${ODIR}/profiling/\n")
            f.write("rm ${IDIR}*\n")

            f.close()
            os.system(f"bash {tmp_fname}")
            os.system(f"rm {tmp_fname}")

class GAL_VOID_Runner(Base_Runner):
    def __init__(self, conf=None):
        super(GAL_VOID_Runner, self).__init__(conf)

    def set_params(self, **kwargs):
        if "input" in self.conf["SHAM"]:
            self.FROMFILE = True
            self.input = self.conf.get("SHAM", "input").strip("\"")
        else:
            self.FROMFILE = False
            self.seedini = self.conf["SHAM"].getint("seedini")

            self.nparams = self.conf["SHAM"].getint("nparams")
            self.param_names = self.conf.get("SHAM", "param_names").strip("\"")
            self.prior_low = list(map(float, self.conf.get("SHAM", "prior_low").split(", ")))
            self.prior_up = list(map(float, self.conf.get("SHAM", "prior_up").split(", ")))

        self.ncats = self.conf["SHAM"].getint("ncats")
        self.nrlzs = self.conf["SHAM"].getint("nrlzs")

        self.ref_num_den = self.conf["SHAM"].getfloat("ref_num_den")
        self.feature = self.conf.get("SHAM", "feature").strip("\"")
        self.z_space = self.conf["SHAM"].getboolean("z_space")

        self.galbase = self.conf.get("SHAM", "outputbase").strip("\"")

        ### Void Finder path
        self.DIVE_exec = self.conf.get("DIVE", "DIVE_EXEC").strip("\"")

        ### load other necessary params
        self.snapdir = self.conf.get("FastPM", "snapdir").strip("\"")
        self.snapbase = self.conf.get("FastPM", "snapbase").strip("\"")
        self.redshifts = list(map(float, self.conf.get("FastPM", "redshifts").split(", ")))

        self.boxsize = self.conf["FastPM"].getfloat("boxsize")
        self.halobase = self.conf.get("ROCKSTAR", "outputbase").strip("\"")

        self.cfgbase = self.conf.get("General", "cfgbase").strip("\"")
        self.cfgsubbase = self.conf.get("General", "cfgsubbase").strip("\"")

        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

    ### TODO
    # def get_sham_params(self):
    #     if self.FROMFILE:
    #         with open(self.input, "r") as f:
    #             line = f.readline().strip("\n")
    #             self.param_names = line.split(" ")[1:]
                
    #         self.param_list = np.loadtxt(self.input, comments="#")
    #         if len(self.param_list) == 1:
    #             self.param_list = self.param_list[np.newaxis, :]
            
    #         self.nparams = self.param_list.shape[1]
    #         if "input" in self.conf["FastPM"]:
    #             cosmo_param_input = self.conf.get("FastPM", "input").strip("\"")
    #             tmp = np.loadtxt()

    def run(self, SHAM_paramfile=None, SHAM_param_list=None, seed1=None, seed2=None, 
            snapname_relic="", haloname="out_0.list"):
        if SHAM_paramfile is not None and os.path.isfile(SHAM_paramfile):
            param_list = np.loadtxt(SHAM_paramfile, comments="#")
        elif SHAM_param_list is not None:
            param_list = SHAM_param_list
        elif self.FROMFILE:
            if os.path.isfile(self.cfgbase+self.input):
                param_list = np.loadtxt(self.cfgbase+self.input)
            else:
                print(f"Cannot find file: {self.cfgbase+self.input}")
                exit()
        else:
            param_rng = np.random.default_rng(seed=seed1)
            param_list = param_rng.uniform(low=self.prior_low, high=self.prior_up, size=(self.ncats, self.nparams))
            f = open("SHAM_list.txt", "w+", encoding="utf-8")
            f.write(f"# {self.param_names}\n")
            np.savetxt(f, param_list)
            f.close()
        
        for zi in (self.redshifts):
            snappath = self.snapdir+self.snapbase+"{}/a_{:.4f}/".format(snapname_relic, 1./(1+zi))
            halo_file = snappath+self.halobase+haloname
            halo = load_rockstar_halo(halo_file, feature=self.feature, z_space=self.z_space)
            for ipar, iparam in enumerate(param_list):
                gal_void_filebase = snappath+self.halobase+self.galbase+f"{ipar}/"
                if not os.path.isdir(gal_void_filebase):
                    os.mkdir(gal_void_filebase)
                for irlz in range(self.nrlzs):
                    ### populate galaxy
                    rlz_rng = np.random.default_rng(seed=seed2+irlz)
                    igal = SHAM_model(iparam, halo, self.feature, self.z_space, self.ref_num_den, rng=rlz_rng)
                    gal_file = gal_void_filebase+f"gal_{self.feature}_rlz{irlz}.txt"
                    fo = open(gal_file, "w+", encoding="utf-8")
                    ### TODO: DIVE does not support header line, maybe change io funs
                    # if z_space:
                    #     fo.write("# x y zrsd\n")
                    # else:
                    #     fo.write("# x y z\n")
                    np.savetxt(fo, igal[:,:3], fmt="%.3f %.3f %.3f")
                    fo.close()
                    ### find voids
                    void_file = gal_void_filebase+f"void_{self.feature}_rlz{irlz}.txt"
                    os.system(f"{self.DIVE_exec} -i {gal_file} -o {void_file} -u {self.boxsize}")

class TPCF_Runner(Base_Runner):
    def __init__(self, conf=None):
        super(TPCF_Runner, self).__init__(conf)
    
    def set_params(self, **kwargs):
        if "input" in self.conf["General"]:
            cfgbase = self.conf.get("General", "cfgbase").strip("\"")
            cosmo_list_file = self.conf.get("General", "input").strip("\"")
            self.ncosmo = len(np.loadtxt(cfgbase+cosmo_list_file))
        else:
            self.ncosmo = self.conf["General"].getint("ncosmo")

        self.pyfcfc_path = self.conf.get("FCFC", "pyFCFC_PATH").strip("\"")
        sys.path.append(self.pyfcfc_path)

        self.cat1 = self.conf.get("FCFC", "cat1").strip("\"")
        self.cat2 = self.conf.get("FCFC", "cat2").strip("\"") 

        if "void" in self.cat1 or "void" in self.cat2:
            ### void size bin info
            self.Rmin = self.conf["FCFC"].getfloat("RVmin")
            self.Rmax = self.conf["FCFC"].getfloat("RVmax")
            self.dRV  = self.conf["FCFC"].getfloat("dRV")

            ### separation bins info
            self.min_sep_inRv = self.conf["FCFC"].getfloat("min_sep_inRv")
            self.max_sep_inRv = self.conf["FCFC"].getfloat("max_sep_inRv")
            self.n_sep_bins  = self.conf["FCFC"].getint("n_sep_bins")
        
        else:
            self.min_sep = self.conf["FCFC"].getfloat("min_sep")
            self.max_sep = self.conf["FCFC"].getfloat("max_sep")
            self.n_sep_bins = self.conf["FCFC"].getint("n_sep_bins")

        if "matter" in self.cat1 or "matter" in self.cat2:
            ### if using approximation to accelerate computing
            if "downsample" in self.conf["FCFC"]:
                self.DOWN_SAMPLE  = True
                self.dsample_rate = self.conf["FCFC"].getfloat("downsample")
            else:
                self.DOWN_SAMPLE  = False

        ### if estimate jackknife error
        if "njk" in self.conf["FCFC"]:
            self.JK = True
            self.njk = self.conf["FCFC"].getint("njk")
        else:
            self.JK = False

        ### output file
        self.outputbase = self.conf.get("FCFC", "outputbase").strip("\"")

        ### input file
        self.nsham_per_cosmo = self.conf["SHAM"].getint("ncats")
        self.nrlzs_per_sham = self.conf["SHAM"].getint("nrlzs")
        self.boxsize = self.conf["FastPM"].getfloat("boxsize")

        self.snapdir = self.conf.get("FastPM", "snapdir").strip("\"")
        self.snapbase = self.conf.get("FastPM", "snapbase").strip("\"")
        self.redshifts = list(map(float, self.conf.get("FastPM", "redshifts").split(", ")))
        self.halobase = self.conf.get("ROCKSTAR", "outputbase").strip("\"")
        self.shambase = self.conf.get("SHAM", "outputbase").strip("\"")
        self.feature = self.conf.get("SHAM", "feature").strip("\"")

        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

    def load_dm_catalog(self, fname):
        tmp = BigFileCatalog(fname, dataset="1/", header="Header")
        data_size = tmp.csize
        dm = tmp['Position'].compute()

        # if apply downsampling
        if self.DOWN_SAMPLE:
            random.seed(0)
            sample_idx = random.sample(list(np.arange(data_size)), int(self.dsample_rate*data_size))
            dm = dm[sample_idx]
        del tmp

        # if estimate jk_err, generate DM random catalog
        if self.JK:
            np.random.seed(0)
            dm_ran = np.random.rand(len(dm)*10, 3)*self.boxsize
            return dm, dm_ran
        else:
            return dm

    def load_void_catalog(self, fname):
        void = np.loadtxt(fname)
        return void
    
    def load_galaxy_catalog(self, fname):
        galaxy = np.loadtxt(fname)
        return galaxy
    
    def get_tpcf_type(self):
        cat1_name = self.cat1
        cat2_name = self.cat2
        tpcf_type = "-".join([cat1_name, cat2_name])
        return tpcf_type
    
    def save_to_file(self, oname, header, data):
        f = open(oname, "w+")
        f.write(header)
        np.savetxt(f, data)
        f.close()
        return None
    
    def get_indices(self, cat, length, nside):
        indices_3d = cat // length
        indices_1d = indices_3d[:,0] + nside*(indices_3d[:,1] + nside*indices_3d[:,2])
        
        return indices_1d

    def declare(self):
        tpcf_type = self.get_tpcf_type()
        print("="*100+"\n")
        print(f"Begin measuring {tpcf_type} two point correlation function ...\n")
        if "matter" in tpcf_type:
            print(f"If downsampling particles: {self.DOWN_SAMPLE}")
            if self.DOWN_SAMPLE:
                print("downsamaple rate: {:.2f}%\n".format(self.dsample_rate*100))
            else:
                print("\n")
        print(f"If estimate Jackknife error: {self.JK}")
        if self.JK:
            print("Jackknife subsamples: {:d}\n".format(self.njk**3))
        else:
            print("\n")
        
        if "void" in tpcf_type:
            print("Void Size bins:\n")
            print("Rmin: {:.2f}\nRmax: {:.2f}\ndR: {:.2f}\n".format(self.Rmin, self.Rmax, self.dRV))
            print("TPCF separation bins:\n")
            print("smin: {:.2f} Rv\nsmax: {:.2f} Rv\nnbins: {:d}\n".format(self.min_sep_inRv, self.max_sep_inRv, self.n_sep_bins))
        else:
            print("TPCF separation bins:\n")
            print("smin: {:.2f} Mpc/h\nsmax: {:.2f} Mpc/h\nnbins: {:d}\n".format(self.min_sep, self.max_sep, self.n_sep_bins))
        print("Catalog info:")
        print("N SHAM per cosmo: {:d}\nN realization per SHAM: {:d}\n".format(self.nsham_per_cosmo, self.nrlzs_per_sham))
        print("="*100)

    def run(self, snapname_relic=""):
        tpcf_type = self.get_tpcf_type()

        if "void" in tpcf_type:
            Rmins = np.arange(self.Rmin, self.Rmax, self.dRV)
            Rmaxs = np.append(Rmins[1:], self.Rmax)

        for zi in self.redshifts:
            if "matter" in tpcf_type:
                # load dm catalog
                snappath = self.snapdir+self.snapbase+"{}/a_{:.4f}/".format(snapname_relic,1./(1.+zi))
                if self.JK:
                    dm, dm_ran = self.load_dm_catalog(snappath)
                else:
                    dm = self.load_dm_catalog(snappath)
            # tmp = BigFileCatalog(snappath, dataset="1/", header="Header")
            # data_size = tmp.csize
            # dm = tmp['Position'].compute()

            # # if apply downsampling
            # if self.DOWN_SAMPLE:
            #     random.seed(0)
            #     sample_idx = random.sample(list(np.arange(data_size)), int(self.dsample_rate*data_size))
            #     dm = dm[sample_idx]
            # del tmp

            # # if estimate jk_err, generate DM random catalog
            # if self.JK:
            #     np.random.seed(0)
            #     dm_ran = np.random.rand(len(dm)*10, 3)*self.boxsize

            for isham in np.arange(self.nsham_per_cosmo):
                outputpath = self.outputbase+self.snapbase+"{}/a_{:.4f}/SHAM{:d}/".format(snapname_relic,1./(1.+zi),isham)
                if not os.path.isdir(outputpath):
                    os.makedirs(outputpath)
                
                for irlz in range(self.nrlzs_per_sham):
                    ### catalog loading part
                    if "void" in tpcf_type:
                        # load void catalog
                        voidpath = snappath+self.halobase+self.shambase+f"{isham}/void_{self.feature}_rlz{irlz}.txt"
                        void = self.load_void_catalog(voidpath)

                    if "galaxy" in tpcf_type:
                        galpath = snappath+self.halobase+self.shambase+f"{isham}/gal_{self.feature}_rlz{irlz}.txt"
                        galaxy = self.load_galaxy_catalog(galpath)
                    
                    ### measurement part
                    if tpcf_type == "matter-galaxy" or tpcf_type == "galaxy-matter":
                        results = get_cross_iso(
                            dm, galaxy, self.boxsize, 
                            bins_info={
                                "min_sep": self.min_sep, 
                                "max_sep": self.max_sep,
                                "nbins": self.n_sep_bins
                            }
                            )
                        
                        oname = outputpath+f"ggl_rlz{irlz}.txt"
                        header = "# sep (Mpc/h) xi_iso\n"
                        data = np.array(results['multipoles'])[0,0,:]
                        self.save_to_file(oname, header, data)

                    if tpcf_type == "matter-void" or tpcf_type == "void-matter":
                        results_list = get_stacked_vprof(
                            dm, void, self.boxsize, Rmins, Rmaxs,
                            min_sep_inRv=self.min_sep_inRv, 
                            max_sep_inRv=self.max_sep_inRv, 
                            nbins=self.n_sep_bins,
                            )

                        if self.JK:
                            # if estimate jk_err, generate void random catalog
                            np.random.seed(0)
                            void_ran = np.random.rand(len(void)*10, 3)*self.boxsize

                            # get subvolume index of each particle/tracer
                            sublength = self.boxsize / self.njk
                            index_dm = self.get_indices(dm, sublength, self.njk)
                            index_dmran = self.get_indices(dm_ran, sublength, self.njk)
                            index_v  = self.get_indices(void, sublength, self.njk)
                            index_vran = self.get_indices(void_ran, sublength, self.njk)
                            
                            # estimate jk err
                            tmp_rls_list = []
                            tmp_xi = []
                            tmp_xi_stacked = []
                            for ivol in range(self.njk*self.njk*self.njk):
                                tmp_rls_list = get_stacked_vprof(dm[~(index_dm == ivol)], void[~(index_v == ivol)], self.boxsize, 
                                                                Rmins, Rmaxs,
                                                                rdm=dm_ran[~(index_dmran == ivol)], rvoid=void_ran[~(index_vran == ivol)], 
                                                                min_sep_inRv=self.min_sep_inRv, 
                                                                max_sep_inRv=self.max_sep_inRv, 
                                                                nbins=self.n_sep_bins
                                                                )

                                tmp_xi.append(tmp_rls_list['xi_iso_list'])
                                tmp_xi_stacked.append(tmp_rls_list['xi_iso_stacked'])

                            tmp_xi = np.asarray(tmp_xi)
                            jk_err_xi = np.sum((tmp_xi - np.mean(tmp_xi, axis=0))**2, axis=0)*((self.njk)**3 - 1)/self.njk**3
                            jk_err_xi = np.sqrt(jk_err_xi)

                            tmp_xi_stacked = np.asarray(tmp_xi_stacked)
                            jk_err_xi_stk = np.sum((tmp_xi_stacked - np.mean(tmp_xi_stacked, axis=0))**2, axis=0)*((self.njk)**3 - 1)/self.njk**3
                            jk_err_xi_stk = np.sqrt(jk_err_xi_stk)

                            for iR, ind_xi_iso in enumerate(results_list["xi_iso_list"]):
                                f = open(outputpath+f"rvbin{iR}_rlz{irlz}.txt", "w+")
                                f.write("# sep (Mpc/h) xi_iso jk_err\n")
                                np.savetxt(f, np.c_[ind_xi_iso, jk_err_xi[iR]])
                                f.close()

                            f = open(outputpath+f"stacked_rlz{irlz}.txt", "w+")
                            f.write("# Rmin={:.2f} Rmax={:.2f} Nbins={:d}\n".format(self.Rmin, self.Rmax, len(Rmins)))
                            f.write("# sep (Rv) xi_iso jk_err\n")
                            np.savetxt(f, np.c_[results_list["xi_iso_stacked"], jk_err_xi_stk])
                            f.close()     

                        else:
                            for iR, ind_xi_iso in enumerate(results_list["xi_iso_list"]):
                                oname = outputpath+f"rvbin{iR}_rlz{irlz}.txt"
                                header = "# sep (Mpc/h) xi_iso\n"
                                data = ind_xi_iso
                                self.save_to_file(oname, header, data)
                                # f = open(outputpath+f"rvbin{iR}_rlz{irlz}.txt", "w+")
                                # f.write("# sep (Mpc/h) xi_iso\n")
                                # np.savetxt(f, ind_xi_iso)
                                # f.close()

                            oname = outputpath+f"stacked_rlz{irlz}.txt"
                            header = "# Rmin={:.2f} Rmax={:.2f} Nbins={:d}\n".format(self.Rmin, self.Rmax, len(Rmins)) + "# sep (Rv) xi_iso\n"
                            data = results_list["xi_iso_stacked"]
                            self.save_to_file(oname, header, data)
                            # f = open(outputpath+f"stacked_rlz{irlz}.txt", "w+")
                            # f.write("# Rmin={:.2f} Rmax={:.2f} Nbins={:d}\n".format(self.Rmin, self.Rmax, len(Rmins)))
                            # f.write("# sep (Rv) xi_iso\n")
                            # np.savetxt(f, results_list["xi_iso_stacked"])
                            # f.close()
                    print(f"Measuring {tpcf_type} 2pcf of {self.snapdir}/{self.snapbase}{snapname_relic}, redshift={zi}, SHAM{isham} Done.", flush=True)


class POWER_Runner(Base_Runner):
    def __init__(self, conf=None):
        super(POWER_Runner, self).__init__(conf)

    def set_params(self, **kwargs):
        ### power spectrum info
        self.Nmesh = self.conf["POWER"].getint("Nmesh") # position to grid
        self.MAS = self.conf.get("POWER", "MAS") # Mass Alignment Scheme

        self.kmin = self.conf["POWER"].getfloat("min_k")
        self.kmax = self.conf["POWER"].getfloat("max_k")
        self.nk   = self.conf["POWER"].getint("nk")
        self.multipoles = list(map(int, self.conf.get("POWER", "multipoles").split(", ")))

        ### output file
        self.outputbase = self.conf.get("POWER", "outputbase").strip("\"")

        ### input file
        if "input" in self.conf["General"]:
            self.cfgbase = self.conf.get("General", "cfgbase").strip("\"")
            self.cosmo_list_file = self.conf.get("General", "input").strip("\"")
            cosmo_file = os.path.join(self.cfgbase, self.cosmo_list_file)
            tmp = np.loadtxt(cosmo_file)
            if len(tmp.shape) == 1:
                self.ncosmo = 1
            else:
                self.ncosmo = len(tmp)
        else:
            self.ncosmo = self.conf["General"].getint("ncosmo")
        self.nsham_per_cosmo = self.conf["SHAM"].getint("ncats")
        self.nrlzs_per_sham = self.conf["SHAM"].getint("nrlzs")
        self.boxsize = self.conf["FastPM"].getfloat("boxsize")

        self.snapdir = self.conf.get("FastPM", "snapdir").strip("\"")
        self.snapbase = self.conf.get("FastPM", "snapbase").strip("\"")
        self.redshifts = list(map(float, self.conf.get("FastPM", "redshifts").split(", ")))
        self.halobase = self.conf.get("ROCKSTAR", "outputbase").strip("\"")
        self.shambase = self.conf.get("SHAM", "outputbase").strip("\"")
        self.feature = self.conf.get("SHAM", "feature").strip("\"")

        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

    def declare(self):
        print("="*100+"\n")
        print("Begin measuring two point correlation function ...\n")
        print("k separation bins:")
        print("kmin: {:.2f} kmax: {:.2f} nbins: {:d}\n".format(self.kmin, self.kmax, self.nk))
        print("Catalog info:")
        print("N SHAM per cosmo: {:d}\nN realization per SHAM: {:d}\n".format(self.nsham_per_cosmo, self.nrlzs_per_sham))
        print("="*100)

    def write_to_file(self, fname, header, data):
        with open(fname, "w+") as f:
            f.write(header)
            if header[-1] != "\n":
                f.write("\n")

            np.savetxt(f, data)
        
        return None

    def run(self, relic="", nthreads=1):
        for zi in self.redshifts:
            snappath = self.snapdir + self.snapbase + "{:s}/a_{:.4f}/".format(relic, 1./(1.+zi))
            for isham in range(self.nsham_per_cosmo):
                outputpath = self.outputbase+self.snapbase+"{:s}/a_{:.4f}/SHAM{:d}/".format(relic,1./(1.+zi),isham)
                if not os.path.isdir(outputpath):
                    os.makedirs(outputpath)

                Pkout_mean = []
                for irlz in range(self.nrlzs_per_sham):
                    galpath = snappath+self.halobase+self.shambase+f"{isham}/gal_{self.feature}_rlz{irlz}.txt"
                    gal = np.loadtxt(galpath)

                    kout, Pkout = get_powerspec_cat(gal, self.Nmesh, self.boxsize, self.MAS, self.kmin, self.kmax, self.nk, self.multipoles, nthreads)
                    fname_irlz = outputpath+f"power_rlz{irlz}.txt"
                    header_irlz = "# k (h/Mpc) " + " ".join([f"Pk{ell}" for ell in range(3)])
                    data_irlz = np.c_[kout, np.asarray(Pkout).T]
                    self.write_to_file(fname_irlz, header_irlz, data_irlz)
                    
                    Pkout_mean.append(Pkout)
                Pkout_mean = np.mean(np.asarray(Pkout_mean), axis=0)

                fname = outputpath+f"power.txt"
                header = "# k (h/Mpc) " + " ".join([f"Pk{ell}" for ell in range(3)])
                data = np.c_[kout, Pkout_mean.T]
                self.write_to_file(fname, header, data)

                # f = open(outputpath+f"power.txt", "w+")
                # f.write("# k (h/Mpc) ")
                # for ell in self.multipoles:
                #     f.write(f"Pk{ell} ")
                # f.write("\n")
                # np.savetxt(f, np.c_[kout, Pkout_mean.T])
                # f.close()