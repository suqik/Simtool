import os

ncosmo = 100
for icosmo in range(ncosmo):
    root = f"/public/share/ace66so15x/suchen/L1000_N1024_validation/cosmo{icosmo}/a_0.7692/rstar/out_0.list"
    if not os.path.isfile(root):
        print(f"{icosmo}")
