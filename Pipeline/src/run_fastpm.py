import os, sys
import configparser

if len(sys.argv) != 3:
    print("Usage: python run_fastpm.py pipeline_conf/fpm_conf nCPUs")
    exit()

conf_file = sys.argv[1]
if not os.path.isfile(os.path.join(os.getcwd(), conf_file)):
    print("Configure file does not exist!")
    exit()

nCPUs = int(sys.argv[2])

if conf_file.split(".")[-1] == "ini":
    conf = configparser.ConfigParser()
    conf.read(conf_file)
    base = conf.get("General","cfgbase").strip("\"")
    subbase = conf.get("General","cfgsubbase").strip("\"")
    ncosmo = conf["General"].getint("ncosmo")

    FastPM_exec = conf.get("FastPM", "FastPM_EXEC").strip("\"")
    cmdbase = f"mpirun -np "+f"{nCPUs} "+FastPM_exec+" "

    os.environ["OMP_NUM_THREADS"] = "1"
    for i in range(ncosmo):
        fpm_cfgpath = base+subbase+f"{i}/fastpm/fpm.lua"
        cmd = cmdbase + fpm_cfgpath
        print(cmd, flush=True)
        os.system(cmd)
