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
gti_dir = "./gti_files"

# Extract spectra for each GTI
spectrum_dir = "./spectra"
os.makedirs(spectrum_dir, exist_ok=True)

for gti_file in sorted(os.listdir(gti_dir)):
    gti_path = os.path.join(gti_dir, gti_file)
    output_spectrum = os.path.join(spectrum_dir, f"spectrum_{os.path.splitext(gti_file)[0]}.fits")
    print(f"Processing {gti_file}...")
    cmd = (
        f"evselect table={table} energycolumn=PI expression='gti({gti_path},TIME)' "
        f"withspectrumset=yes spectrumset={output_spectrum} "
        f"specchannelmin=0 specchannelmax=20479 spectralbinsize=5"
    )
    subprocess.run(cmd, shell=True, check=True)

print(f"Spectra saved in {spectrum_dir}.")