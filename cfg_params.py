### This file contains parameter dictionarys for different parameter files. ###
param_list = ['fastpm', 'preproc', 'gluons']
### FastPM
fastpm_params={'nc':                       256,
                'boxsize':                 400.0,
                'time_step':               "linspace(0.01, 1, 40)",
                'output_redshifts':        "{0.0}",
                'Omega_m':                 0.3156,
                'hubble':                  0.6727,
                'read_powerspectrum':      "\"powerspec.txt\"",
                'linear_density_redshift': 0.0,
                'random_seed':             100,
                'force_mode':              "\"fastpm\"",
                'pm_nc_factor':            2,
                'np_alloc_factor':         3,
                'write_snapshot':          "\"fastpm\"",
                'write_powerspectrum':     None
}
fpm_int_keys = ['nc', 'random_seed', 'pm_nc_factor', 'np_alloc_factor']
fpm_float_keys = ['boxsize', 'Omega_m', 'hubble', 'linear_density_redshift']
fpm_special_keys = ['time_step', 'output_redshifts']

### LightCone Generator
preproc_params={'SIMPATH': "/home/suqikuai777/Program/SBI_Void/simulator/fastpm/catalog/match_Nbody",
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
extra_params = ['ofile', 'bsfile', 's', 'b', 'w']
lc_int_keys = ['SIMNUM', 'NPART', 'FILENUM', 'NMEMBERS', 'LIST']
lc_float_keys = ['LINKLENGTH', 'OVERLAP']

### Raytracing 
gluons_params = {'OMEGAX': 0.6727,
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
ray_int_keys = ['NUMSLC', 'NUMDEC', 'NUMRA', 'NGRID']
ray_float_keys = ['OMEGAX', 'OMEGAM', 'HUBBLE', 'RNGDEC', 'RNGRA', 'CTRDEC', 'CTRRA', 'SMOOTH', 'ZSOURC']