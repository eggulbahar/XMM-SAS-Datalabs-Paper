from pysas.wrapper import Wrapper as w
import os.path
from os import path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table



wdir=os.getcwd()
home = os.path.expanduser('~')
os.environ['SAS_CCFPATH'] = f'{home}/data/user/pub'
inargs = [f'sas_ccf={wdir}/ccf.cif', f'sas_odf={wdir}/3553_0841890201_SCX00000SUM.SAS', f'workdir={wdir}']
w('startsas', inargs).run()

table = wdir + "/PN_clean_evt.fits"

with fits.open(table) as hdul:
    header = hdul[1].header 
    obs_start = header.get("TSTART")
    obs_end = header.get("TSTOP")
    print(f"TIME-OBS: {obs_start}")
    print(f"TIME-END: {obs_end}")


# Create GTI files for each pulse period
gti_dir = "./gti_files"
os.makedirs(gti_dir, exist_ok=True)

print("Creating 283.44-second GTI files...")
start_time = obs_start

while start_time < obs_end:
    stop_time = start_time + 283.44
    gti_file = os.path.join(gti_dir, f"gti_{start_time:.3f}_{stop_time:.3f}.fits")
    cmd = f"tabgtigen table={table} expression='TIME >= {start_time} && TIME < {stop_time}' gtiset={gti_file}"
    subprocess.run(cmd, shell=True, check=True)
    start_time = stop_time

print(f"GTI files created in {gti_dir}.")