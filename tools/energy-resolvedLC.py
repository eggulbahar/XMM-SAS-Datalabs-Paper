from pysas.wrapper import Wrapper as w
import os.path
from os import path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter  # For formatting axis labels


wdir=os.getcwd()
home = os.path.expanduser('~')
os.environ['SAS_CCFPATH'] = f'{home}/data/user/pub'
inargs = [f'sas_ccf={wdir}/ccf.cif', f'sas_odf={wdir}/3553_0841890201_SCX00000SUM.SAS', f'workdir={wdir}']
w('startsas', inargs).run()

table = wdir + "/PN_clean_evt.fits"

# Avoiding pile-up regions
rawX1src= 32
rawX2src = 44
rawX3src = 36
rawX4src = 40

rawX1bkg = 3 #RAWX1 of the bkg rectangle region (in pixels)
rawX2bkg = 5 #RAWX2 of the bkg rectangle region (in pixels)

energy_ranges= [500, 3000, 6000, 8000, 10000]

EresolvedLC=[]

for i in range(len(energy_ranges) - 1):
    
    # Define E ranges
    e_min = energy_ranges[i]  # Start of the E range 
    e_max = energy_ranges[i+1]  # End of the E range
    
    # Extract the source region lightcurve

    # Define some parameters for filtering the event file and define the lightcurve binning

    q_flag       = "#XMMEA_EP" # Quality flag for EPIC pn
    n_pattern    = 4           # Pattern selection
    pn_pi_min    = e_min        # Low energy range eV
    pn_pi_max    = e_max     # High energy range eV
    lc_bin       = 283         # Lightcurve bin in secs

    # Define the output ligthcurve file name

    in_LCSRCFile = wdir+'/PN_source_lightcurve_raw_'+str(e_min)+'to'+str(e_max)+'eV_bin283sec.lc'   # Name of the output source lightcurve

    # SAS Command
    cmd        = "evselect" # SAS task to be executed                  

    # Arguments of SAS Command
    expression = f'{q_flag}&&(PATTERN<={n_pattern})&&(RAWX in [{rawX1src}:{rawX3src}] || RAWX in [{rawX4src}:{rawX2src}])&&(PI in [{pn_pi_min}:{pn_pi_max}])'  # event filter expression
    inargs     = [f'table={table}','energycolumn=PI','withrateset=yes',f'rateset={in_LCSRCFile}',
                  f'timebinsize={lc_bin}','maketimecolumn=yes','makeratecolumn=yes',f'expression={expression}']


    w(cmd, inargs).run()
    
    # Extract the background region lightcurve

    q_flag       = "#XMMEA_EP" # Quality flag for EPIC pn
    n_pattern    = 4           # Pattern selection
    pn_pi_min    = e_min        # Low energy range eV
    pn_pi_max    = e_max      # High energy range eV
    lc_bin       = 283         # Lightcurve bin in secs

    # Define the output ligthcurve file name

    in_LCBKGFile = wdir+'/PN_lightcurve_background_raw_'+str(e_min)+'to'+str(e_max)+'eV_bin283sec.lc'   # Name of the output source lightcurve

    cmd        = "evselect" # SAS task to be executed                  

    # Arguments of SAS Command
    expression = f'{q_flag}&&(PATTERN<={n_pattern})&&(RAWX in [{rawX1bkg}:{rawX2bkg}])&&(PI in [{pn_pi_min}:{pn_pi_max}])'  # event filter expression
    inargs     = [f'table={table}','energycolumn=PI','withrateset=yes',f'rateset={in_LCBKGFile}',
                  f'timebinsize={lc_bin}','maketimecolumn=yes','makeratecolumn=yes',f'expression={expression}']
    w(cmd, inargs).run()
    
    # Extract the corrected lightcurve

    in_LCFile = wdir+'/PN_lccorr_'+str(e_min)+'to'+str(e_max)+'eV_bin283sec.lc'   # Name of the output corrected lightcurve
    cmd        = "epiclccorr" # SAS task to be executed                  

    # Arguments of SAS Command
    inargs     = [f'eventlist={table}',f'srctslist={in_LCSRCFile}',f'outset={in_LCFile}',
                  f'bkgtslist={in_LCBKGFile}','withbkgset=yes','applyabsolutecorrections=yes']
    w(cmd, inargs).run()
    
    EresolvedLC.append(in_LCFile)

print(f'All the energy-resolved light curves produced: {EresolvedLC}')