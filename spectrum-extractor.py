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



T_obs=[58606.95, 58607.6, 58607.78, 58608.2]
tt_times=[]
for i, t in enumerate(T_obs):
    relative_seconds =(t-58607)* 86400
    time_difference = (58607 - 50814) * 86400
    tt_times.append(relative_seconds+time_difference)
print(tt_times)

wdir=os.getcwd()
home = os.path.expanduser('~')
os.environ['SAS_CCFPATH'] = f'{home}/data/user/pub'
inargs = [f'sas_ccf={wdir}/ccf.cif', f'sas_odf={wdir}/3553_0841890201_SCX00000SUM.SAS', f'workdir={wdir}']
w('startsas', inargs).run()

table = wdir + "/PN_clean_evt.fits"

grouped_spectra=[]

for i in range(len(tt_times) - 1):
    
    # Define time ranges
    time_min = tt_times[i]  # Start of the time range (in seconds, for example)
    time_max = tt_times[i+1]  # End of the time range
    rawX1src= 32
    rawX2src = 44
    rawX3src = 36
    rawX4src = 40
    
    spectrumset= wdir+ '/PN_source_spectrum_raw_'+str(tt_times[i])+'_'+str(rawX1src)+'-'+str(rawX3src)+'_'+str(rawX4src)+'-'+str(rawX2src)+'.fits'

    # SAS Command
    cmd        = "evselect" # SAS task to be executed                  

    
    expression = f'(FLAG==0) && (PATTERN<=4) && (RAWX in [{rawX1src}:{rawX3src}] || RAWX in [{rawX4src}:{rawX2src}]) && (TIME >= {time_min}) && (TIME <= {time_max})'

    # Arguments of SAS Command
    inargs = [f'table={table}', 'withspectrumset=yes', f'spectrumset={spectrumset}',
              'energycolumn=PI', 'spectralbinsize=5', 'withspecranges=yes',
              'specchannelmin=0', 'specchannelmax=20479', f'expression={expression}']
    w(cmd, inargs).run()
    
    
    cmd        = "backscale" # SAS task to be executed                  
    inargs     = [f'spectrumset={spectrumset}',f'badpixlocation={table}']
    w(cmd, inargs).run()
    
    in_RESPFile = wdir+'/PN_'+str(tt_times[i])+'_'+str(rawX1src)+'-'+str(rawX3src)+'_'+str(rawX4src)+'-'+str(rawX2src)+'.rmf' 
    cmd        = "rmfgen" # SAS task to be executed                  
    inargs     = [f'spectrumset={spectrumset}',f'rmfset={in_RESPFile}']
    w(cmd, inargs).run()
    
    in_ARFFile = wdir+'/PN_'+str(tt_times[i])+'_'+str(rawX1src)+'-'+str(rawX3src)+'_'+str(rawX4src)+'-'+str(rawX2src)+'.arf'
    cmd        = "arfgen" # SAS task to be executed                  

    inargs     = [f'spectrumset={spectrumset}',f'arfset={in_ARFFile}',
              'withrmfset=yes',f'rmfset={in_RESPFile}',f'badpixlocation={table}','detmaptype=psf', 'applyabsfluxcorr=yes']
    w(cmd, inargs).run()
    
    
    in_GRPFile = wdir+'/PN_spectrum_grp_'+str(tt_times[i])+'_'+str(rawX1src)+'-'+str(rawX3src)+'_'+str(rawX4src)+'-'+str(rawX2src)+'.fits'  
    cmd        = "specgroup" # SAS task to be executed                  

    inargs     = [f'spectrumset={spectrumset}','mincounts=25','oversample=3',
              f'rmfset={in_RESPFile}',f'arfset={in_ARFFile}',
              f'groupedset={in_GRPFile}']
    w(cmd, inargs).run()
    
    grouped_spectra.append(in_GRPFile)
    

print(grouped_spectra)