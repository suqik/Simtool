import sys
sys.path.append("/home/suqikuai777/Program/SBI_Void/simulator/fastpm")
from Simtool.cfg_params import *
import numpy as np
import warnings

Usage_message = "Usage: python make-cfgs.py cfg_type [options].\n\
cfg_type (case insensitive): fastpm  preproc  gluons \n\
\n\
Options for fastpm cfg: \n\
    -nc                        number of particle \n\
    -boxsize                   length of box in Mpc/h \n\
    -time_step                 array of time step (in scale factor) of simulation \n\
                                example: \"linspace(0.01, 1, 40)\" \n\
    -output_redshifts          dictionary of output redshifts. \n\
                                example: \"{0.0}\" \n\
    -Omega_m                   total matter density parameter \n\
    -hubble                    Hubble parameter in km/s/Mpc \n\
    -read_powerspectrum        path of input power spectrum \n\
    -read_lineark              path of input initial matter distribution \n\
    -linear_density_redshift   redshift of input power spectrum \n\
    -random_seed               random seed\n\
    -force_mode                Force mode used in FastPM \n\
    -pm_nc_factor              parameter B in FastPM \n\
    -np_alloc_factor           ratio of allocated memory for each particle to its fiducial \n\
    -write_snapshot            path of output snapshots. None if do not output \n\
    -write_powerspectrum       path of output power spectrum at each time step. None if do not output\n\
    -ofile                     path of configuration file. Default is the work dictionary, named `test.lua`"

cfg_type = sys.argv[1]

if cfg_type.lower() not in param_list:
    print(Usage_message)
    if cfg_type != '-h' and cfg_type != '--help':
        print(f"First param needs to be type of parameter file. Cannot recognize type of parameter file: {cfg_type}.")
    exit(1)

if cfg_type.lower() == 'fastpm':
    ofile = "test.lua"
    curr_params = fastpm_params.copy()

    for idx, param in enumerate(sys.argv):
        if param[0] == '-':
            value = sys.argv[idx+1]
            if param[1:] not in curr_params.keys() and param != '-ofile':
                warnings.warn(f"Cannot recognize key {param}. Ignored.")
            else:
                param = param[1:]
                if param == 'ofile':
                    ofile = value
                elif param in fpm_int_keys:
                    curr_params[param] = int(value)
                elif param in fpm_float_keys:
                    curr_params[param] = float(value)
                elif param in fpm_special_keys:
                    curr_params[param] = value                
                else:
                    curr_params[param] = "\""+(value)+"\""

    f = open(ofile, 'w+', encoding='utf-8')
    for key in curr_params.keys():
        if curr_params[key] is not None:
            if key == 'hubble':
                f.write(f"h = {curr_params[key]}\n")
            else:
                f.write(f"{key} = {curr_params[key]}\n")

    f.close()

if cfg_type.lower() == 'preproc':
    ofile = "preproc.param"
    curr_params = preproc_params.copy()
    if '-bsfile' not in sys.argv:
        if not ('-s' in sys.argv and '-b' in sys.argv and '-w' in sys.argv):
            print("Need box-snap list file as input, or parameters to generate the file.")
            exit(1)
        else:
            bsf = False
    else:
        bsf = True

    for idx, param in enumerate(sys.argv):
        if param[0] == '-':
            value = sys.argv[idx+1]
            if param[1:] not in curr_params.keys() and param[1:] not in extra_params:
                warnings.warn(f"Cannot recognize key {param}. Ignored.")
            else:
                param = param[1:]
                if param == 'ofile':
                    ofile = value
                elif param == 'bsfile':
                    if '-s' in sys.argv or '-b' in sys.argv or '-w' in sys.argv:
                        warnings.warn('Box-snap list file is priority.')
                    else:
                        bsfile = value
                elif not bsf and param == 's':
                    nslice = int(value)
                elif not bsf and param == 'b':
                    nbox = int(value)
                elif not bsf and param == 'w':
                    width = float(value)
                elif param in curr_params.keys():
                    if param in lc_int_keys:
                        curr_params[param] = int(value)
                    elif param in lc_float_keys:
                        curr_params[param] = float(value)
                    else:
                        curr_params[param] = (value)

    f = open(ofile, 'w+', encoding='utf-8')
    for key in curr_params.keys():
        if curr_params[key] is not None:
            f.write(f"{key} {curr_params[key]}\n")
    if bsf:
        bslist = np.loadtxt(bsfile, comments='#')
        np.savetxt(f, bslist, fmt="%d %d %.2f")
    else:
        buff = []
        for i in range(nslice):
            if i != 0 and i//nbox != (i-1)//nbox:
                f.writelines(buff[::-1])
                buff = []
            buff.append("{} {} {:.2f}\n".format(i//nbox, i%nbox, width))
        if len(buff) != 0:
            f.writelines(buff[::-1])
        # for i in range(nslice):
        #     f.write("{} {} {:.2f}\n".format(i//nbox, nbox-1-i%nbox, width))

    f.close()

if cfg_type.lower() == 'gluons':
    ofile = "raytrace.param"
    curr_params = gluons_params.copy()

    for idx, param in enumerate(sys.argv):
        if param[0] == '-':
            value = sys.argv[idx+1]
            if param[1:] not in curr_params.keys() and param[1:] not in extra_params:
                warnings.warn(f"Cannot recognize key {param}. Ignored.")
            else:
                param = param[1:]
                if param == 'ofile':
                    ofile = value
                elif param in curr_params.keys():
                    if param in ray_int_keys:
                        curr_params[param] = int(value)
                    elif param in ray_float_keys:
                        curr_params[param] = float(value)
                    else:
                        curr_params[param] = (value)

    f = open(ofile, 'w+', encoding='utf-8')
    for key in curr_params.keys():
        if curr_params[key] is not None:
            f.write(f"{key} {curr_params[key]}\n")

    f.close()

if cfg_type.lower() == 'fcfc':
    ofile = "fcfc_2pt_box.conf"
    curr_params = fcfc_params.copy()

    for idx, param in enumerate(sys.argv):
        if param[0] == '-':
            value = sys.argv[idx+1]
            if param[1:] not in curr_params.keys() and param[1:] not in extra_params:
                warnings.warn(f"Cannot recognize key {param}. Ignored.")
            else:
                param = param[1:]
                if param == 'ofile':
                    ofile = value
                elif param in curr_params.keys():
                    if param in fcfc_int_keys:
                        curr_params[param] = int(value)
                    elif param in fcfc_float_keys:
                        curr_params[param] = float(value)
                    else:
                        curr_params[param] = (value)

    f = open(ofile, 'w+', encoding='utf-8')
    for key in curr_params.keys():
        if curr_params[key] is not None:
            f.write(f"{key} = {curr_params[key]}\n")

    f.close()    