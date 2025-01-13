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

# Avoiding pile-up regions
rawX1src= 32
rawX2src = 44
rawX3src = 36
rawX4src = 40

for gti_file in sorted(os.listdir(gti_dir)):
    gti_path = os.path.join(gti_dir, gti_file)
    output_spectrum = os.path.join(spectrum_dir, f"spectrum_{os.path.splitext(gti_file)[0]}.fits")
    print(f"Processing {gti_file}...")
    cmd = "evselect"
    expression = f'(FLAG==0) && (PATTERN<=4) && (RAWX in [{rawX1src}:{rawX3src}] || RAWX in [{rawX4src}:{rawX2src}]) && (gti({gti_path},TIME))' 
    inargs = [f'table={table}', 'withspectrumset=yes', f'spectrumset={output_spectrum}', 'energycolumn=PI', 'spectralbinsize=5', 'withspecranges=yes', 'specchannelmin=0', 'specchannelmax=20479', f'expression={expression}']
    w(cmd, inargs).run()

print(f"Spectra saved in {spectrum_dir}.")