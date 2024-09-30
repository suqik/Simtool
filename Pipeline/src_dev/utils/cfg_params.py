### This file contains parameter dictionarys for different parameter files. ###
param_list = ['fastpm', 'rockstar', 'preproc', 'gluons', 'fcfc']
### FastPM
fastpm_default={'nc':                       256,
                'boxsize':                 400.0,
                'time_step':               "linspace(0.01, 1, 40)",
                'output_redshifts':        "{0.0}",
                'Omega_m':                 0.3156,
                'hubble':                  0.6727,
                'read_lineark':            None,
                'read_powerspectrum':      "\"powerspec.txt\"",
                'remove_cosmic_variance':  False, 
                'inverted_ic':             False,
                'linear_density_redshift': 0.0,
                'random_seed':             100,
                'force_mode':              "\"fastpm\"",
                'pm_nc_factor':            2,
                'np_alloc_factor':         3,
                'write_snapshot':          None,
                'write_powerspectrum':     None,
                'write_rfof':              None,
                'rfof_nmin':               None
}
# fpm_int_keys = ['nc', 'random_seed', 'pm_nc_factor', 'np_alloc_factor', 'rfof_nmin']
# fpm_float_keys = ['boxsize', 'Omega_m', 'hubble', 'linear_density_redshift']
# fpm_special_keys = ['time_step', 'output_redshifts']

### Rockstar
rockstar_default={'FILE_FORMAT': repr("GADGET2"),
                 'GADGET_LENGTH_CONVERSION': 1,
                 'GADGET_MASS_CONVERSION': "1e+10",
                 'FORCE_RES': 0.003,
                 'PARALLEL_IO': 1,
                 'INBASE': None,
                 'FILENAME': repr("snapshot_000.<block>"),
                 'NUM_BLOCKS': 64,
                 'NUM_READERS': 32,
                 'NUM_WRITERS': 32,
                 'FORK_READERS_FROM_WRITERS': 1,
                 'FORK_PROCESSORS_PER_MACHINE': 32,
                 'OUTBASE': None}

### LightCone Generator
preproc_default={'SIMPATH': "/home/suqikuai777/Program/SBI_Void/simulator/fastpm/catalog/match_Nbody",
                'SIMBASE': "box",
                'SNAPBASE': "gadget/snapshot",
                'SIMNUM': 3,
                'NPART':  256,
                'FILENUM': 1,
                'LINKLENGTH': 0.18,
                'NMEMBERS': 35,
                'OVERLAP': 3000.0,
                'SLICEPATH': "/home/suqikuai777/Program/SBI_Void/simulator/fastpm/catalog/match_Nbody/LC",
                'HALOPATH': "/home/suqikuai777/Program/SBI_Void/simulatot/fastpm/catalog/match_Nbody/HALO",
                'LIST': 1
}
# extra_params = ['ofile', 'bsfile', 's', 'b', 'w']
# lc_int_keys = ['SIMNUM', 'NPART', 'FILENUM', 'NMEMBERS', 'LIST']
# lc_float_keys = ['LINKLENGTH', 'OVERLAP']

### Raytracing 
gluons_default = {'OMEGAX': 0.6727,
                 'OMEGAM': 0.3156,
                 'HUBBLE': 0.6727,
                 'NUMSLC': 23,
                 'PTHSLC': "/home/suqikuai777/Program/SBI_Void/simulator/fastpm/catalog/match_Nbody/LC/slice",
                 'PTHRAY': "NULL",
                 'OUTPUT': "/home/suqikuai777/Program/SBI_Void/simulator/fastpm/catalog/match_Nbody/lens/fov_10sq",
                 'RNGDEC': 36000.0,
                 'RNGRA':  36000.0,
                 'NUMDEC': 512,
                 'NUMRA':  512,
                 'CTRDEC': 0.5,
                 'CTRRA':  0.5,
                 'NGRID':  16384,
                 'SMOOTH': 30.0,
                 'ZSOURC': 0.75
}
# ray_int_keys = ['NUMSLC', 'NUMDEC', 'NUMRA', 'NGRID']
# ray_float_keys = ['OMEGAX', 'OMEGAM', 'HUBBLE', 'RNGDEC', 'RNGRA', 'CTRDEC', 'CTRRA', 'SMOOTH', 'ZSOURC']

### 2PCF
fcfc_default = {'CATALOG': "[simu.cat, rand.cat]",
               'CATALOG_LABEL': "[D, R]",
               'ASCII_FORMATTER': "[\"%f %f %f\", \"%f %f %f\"]",
               'POSITION': "[$1, $2, $3, $1, $2, $3]",
               'SELECTION': None,
               'BOX_SIZE': 400.0,
               'BINNING_SCHEME': 0,
               'PAIR_COUNT': "[DD, DR, RR]",
               'PAIR_COUNT_FILE': "[pair.dd, pair.dr, pair.rr]",
               'CF_ESTIMATOR': "(DD - 2DR + RR)/RR",
               'CF_OUTPUT_FILE': "2pcf.dat",
               'SEP_BIN_FILE': None,
               'SEP_BIN_MIN': 0.0,
               'SEP_BIN_MAX': 100.0,
               'SEP_BIN_SIZE': 1.0,
               'OUTPUT_FORMAT': 1,
               'OVERWRITE': 1
}
fcfc_int_keys = ['BINNING_SCHEME', 'OUTPUT_FORMAT', 'OVERWRITE']
fcfc_float_keys = ['BOX_SIZE', 'SEP_BIN_MIN', 'SEP_BIN_MAX', 'SEP_BIN_SIZE']