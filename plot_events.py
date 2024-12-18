"""
Author: Esin G. Gulbahar
Last update: 13/12/2024

This Python script provides functions to analyze and visualize X-ray astronomy event data from FITS files. It includes the following functionalities:

1. plot_events: 
   - Visualizes event data from one or more FITS event files.
   - Supports applying filters based on event properties (e.g., pattern, PI range, flags, CCD number).
   - Offers options for plotting raw and filtered events side-by-side.
   - Supports multiple coordinate systems (e.g., detector or sky coordinates).

2. gti_filtered_image:
   - Generates side-by-side comparisons of raw and Good Time Interval (GTI) filtered event images.
   - Facilitates visualization of the effects of GTI filtering on event data.
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from matplotlib.colors import LogNorm
import re

# Function to plot events from one or more event files
def plot_events(eventfile, coordinates, plot_both=False, mask=False, expression= None, pattern=None, flag=None, pi_min=None, pi_max=None, ccdnr=None):
        
    # Check if both raw and filtered events should be plotted together
    if mask:
        if plot_both:
            evts=len(eventfile)
            plt.figure(figsize=(20,7))
            
            pl=1 # Initialize plot counter
           
            hdu_list = fits.open(eventfile, memmap=True)
            evt_data = Table(hdu_list[1].data)
            
            
            if expression is not None:

                # Define regular expressions to match the pattern and values
                pattern_regex = r'PATTERN\s*<=\s*(\d+)'
                pi_range_regex = r'PI\s+in\s+\[([\d.]+):([\d.]+)\]'
                flag_regex = r'FLAG\s*==\s*(\d+)'
                ccd_regex = r'CCDNR\s*==\s*(\d+)'

                # Find matches in the expression
                pattern_match = re.search(pattern_regex, expression)
                pi_range_match = re.search(pi_range_regex, expression)
                flag_match = re.search(flag_regex, expression)
                ccd_match = re.search(ccd_regex, expression)

                # Extract values from matches
                pattern = int(pattern_match.group(1)) if pattern_match else 4
                pi_min = float(pi_range_match.group(1)) if pi_range_match else None
                pi_max = float(pi_range_match.group(2)) if pi_range_match else None
                flag = int(flag_match.group(1)) if flag_match else 0
                ccdnr = int(ccd_match.group(1)) if ccd_match else None
                
            if ccdnr is not None:
                mask = ((evt_data['PATTERN'] <= pattern) &
                    (evt_data['FLAG'] == flag) &
                    (evt_data['PI'] >= pi_min) &
                    (evt_data['PI'] <= pi_max) &
                    (evt_data['CCDNR'] == ccdnr))
                    
                ccdnrmask = ((evt_data['CCDNR'] == ccdnr))
            else:
        
                # Apply filters to the data
                mask = ((evt_data['PATTERN'] <= pattern) &
                        (evt_data['FLAG'] == flag) &
                        (evt_data['PI'] >= pi_min) &
                        (evt_data['PI'] <= pi_max))
            
            # Print information about events in the event file and filtered event file
            print("Events in event file" + " " + eventfile + ": " + str(len(evt_data)) + "\n")
            print("Events in filtered event file" + " " + eventfile + ": " + str(np.sum(mask)) + "\n")

            # Create Events image    

            if coordinates == 'DET':
                x_data = evt_data['DETX']
                y_data = evt_data['DETY']
                x_label = 'DETX'
                y_label = 'DETY'
            else:
                x_data = evt_data['X']
                y_data = evt_data['Y']
                x_label = 'X'
                y_label = 'Y'
                
            if ccdnr is not None:
                xmax=np.amax(x_data[ccdnrmask])
                xmin=np.amin(x_data[ccdnrmask])
                xmid=abs(xmax-xmin)/2.+xmin
                ymax=np.amax(y_data[ccdnrmask])
                ymin=np.amin(y_data[ccdnrmask])
                
                xbin_size=80
                ybin_size=80
                NBINS = (int(abs(xmax-xmin)/xbin_size),int(abs(ymax-ymin)/ybin_size))

                plt.subplot(1, 2, pl)

                img_zero_mpl = plt.hist2d(x_data[ccdnrmask], y_data[ccdnrmask], NBINS, cmap='GnBu', norm=LogNorm())
            else:
                xmax=np.amax(x_data)
                xmin=np.amin(x_data)
                xmid=abs(xmax-xmin)/2.+xmin
                ymax=np.amax(y_data)
                ymin=np.amin(y_data)
            
                xbin_size=80
                ybin_size=80
                NBINS = (int(abs(xmax-xmin)/xbin_size),int(abs(ymax-ymin)/ybin_size))

                plt.subplot(1, 2, pl)

                img_zero_mpl = plt.hist2d(x_data, y_data, NBINS, cmap='GnBu', norm=LogNorm())

            cbar = plt.colorbar(ticks=[10.,100.,1000.])
            cbar.ax.set_yticklabels(['10','100','1000'])

            plt.title(eventfile)
            plt.xlabel(x_label)
            plt.ylabel(y_label)

            pl=pl+1  # Increment plot counter

            # Create Filtered Events image

            xmax=np.amax(x_data[mask])
            xmin=np.amin(x_data[mask])
            xmid=abs(xmax-xmin)/2.+xmin
            ymax=np.amax(y_data[mask])
            ymin=np.amin(y_data[mask])
            xbin_size=80
            ybin_size=80
            NBINS = (int(abs(xmax-xmin)/xbin_size),int(abs(ymax-ymin)/ybin_size))

            plt.subplot(1, 2, pl)

            img_zero_mpl = plt.hist2d(x_data[mask], y_data[mask], NBINS, cmap='GnBu', norm=LogNorm())

            cbar = plt.colorbar(ticks=[10.,100.,1000.])
            cbar.ax.set_yticklabels(['10','100','1000'])

            plt.title(eventfile)
            plt.xlabel(x_label)
            plt.ylabel(y_label)

            # Add text with filter criteria
            if ccdnr is not None:
                txt=("PATTERN <= " + str(pattern) +
                        " : " + str(pi_min) + " <= E(eV) <= " + str(pi_max) +
                        " : " + " FLAG = " + str(flag) +
                        " : " + " CCDNR = " + str(ccdnr))
            else:
                txt=("PATTERN <= " + str(pattern) +
                        " : " + str(pi_min) + " <= E(eV) <= " + str(pi_max) +
                        " : " + " FLAG = " + str(flag))
                
            plt.text(xmid, ymin+0.1*abs(ymax-ymin), txt, ha='center')

            pl=pl+1

            hdu_list.close()


        else:
            plt.figure(figsize=(10,7))
            hdu_list = fits.open(eventfile, memmap=True)
            evt_data = Table(hdu_list[1].data)
            prihdu   = hdu_list[1].header

            # Read some information from keywords to be used later on

            if ('INSTRUME' in prihdu):
                ins = prihdu['INSTRUME']
                print("Looking into instrument: "+ins+" \n")
            if ('EXPIDSTR' in prihdu):
                expid = prihdu['EXPIDSTR']
                print("Looking at exposure: "+expid+" \n")

            # Check number of event in initial event file

            print("Events in event file" + " " + eventfile + ": " + str(len(evt_data)) + "\n")
            
            if expression is not None:
                # Define regular expressions to match the pattern and values
                pattern_regex = r'PATTERN\s*<=\s*(\d+)'
                pi_range_regex = r'PI\s+in\s+\[([\d.]+):([\d.]+)\]'
                flag_regex = r'FLAG\s*==\s*(\d+)'
                ccd_regex = r'CCDNR\s*==\s*(\d+)'

                # Find matches in the expression
                pattern_match = re.search(pattern_regex, expression)
                pi_range_match = re.search(pi_range_regex, expression)
                flag_match = re.search(flag_regex, expression)
                ccd_match = re.search(ccd_regex, expression)

                # Extract values from matches
                pattern = int(pattern_match.group(1)) if pattern_match else 4
                pi_min = float(pi_range_match.group(1)) if pi_range_match else None
                pi_max = float(pi_range_match.group(2)) if pi_range_match else None
                flag = int(flag_match.group(1)) if flag_match else 0
                ccdnr = int(ccd_match.group(1)) if ccd_match else None

            # Apply filters to the data
            if ccdnr is not None:
                mask = ((evt_data['PATTERN'] <= pattern) &
                        (evt_data['FLAG'] == flag) &
                        (evt_data['PI'] >= pi_min) &
                        (evt_data['PI'] <= pi_max) &
                        (evt_data['CCDNR'] == ccdnr))
            else:
                mask = ((evt_data['PATTERN'] <= pattern) &
                        (evt_data['FLAG'] == flag) &
                        (evt_data['PI'] >= pi_min) &
                        (evt_data['PI'] <= pi_max))

            if coordinates == 'DET':
                x_data = evt_data['DETX']
                y_data = evt_data['DETY']
                x_label = 'DETX'
                y_label = 'DETY'
            else:
                x_data = evt_data['X']
                y_data = evt_data['Y']
                x_label = 'X'
                y_label = 'Y'

            xmax=np.amax(x_data[mask])
            xmin=np.amin(x_data[mask])
            xmid=abs(xmax-xmin)/2.+xmin
            ymax=np.amax(y_data[mask])
            ymin=np.amin(y_data[mask])
            xbin_size=80
            ybin_size=80
            NBINS = (int(abs(xmax-xmin)/xbin_size),int(abs(ymax-ymin)/ybin_size))

            # Plot image


            img_zero_mpl = plt.hist2d(x_data[mask], y_data[mask], NBINS, cmap='GnBu', norm=LogNorm())

            cbar = plt.colorbar(ticks=[10.,100.,1000.])
            cbar.ax.set_yticklabels(['10','100','1000'])

            plt.title(eventfile)
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            if ccdnr is not None:
                txt=("PATTERN <= " + str(pattern) +
                            " : " + str(pi_min) + " <= E(eV) <= " + str(pi_max) +
                            " : " + " FLAG = " + str(flag) +
                            " : " + " CCDNR = " + str(ccdnr))
            else:
                txt=("PATTERN <= " + str(pattern) +
                            " : " + str(pi_min) + " <= E(eV) <= " + str(pi_max) +
                            " : " + " FLAG = " + str(flag))
                
            plt.text(xmid, ymin+0.1*abs(ymax-ymin), txt, ha='center')
            plt.show()
    
    else:
        plt.figure(figsize=(10,7))
        hdu_list = fits.open(eventfile, memmap=True)
        evt_data = Table(hdu_list[1].data)

        if coordinates == 'DET':
                x_data = evt_data['DETX']
                y_data = evt_data['DETY']
                x_label = 'DETX'
                y_label = 'DETY'
        else:
                x_data = evt_data['X']
                y_data = evt_data['Y']
                x_label = 'X'
                y_label = 'Y'

        xmax=np.amax(x_data)
        xmin=np.amin(x_data)
        xmid=abs(xmax-xmin)/2.+xmin
        ymax=np.amax(y_data)
        ymin=np.amin(y_data)
        xbin_size=80
        ybin_size=80
        NBINS = (int(abs(xmax-xmin)/xbin_size),int(abs(ymax-ymin)/ybin_size))

        img_zero_mpl = plt.hist2d(x_data, y_data, NBINS, cmap='GnBu', norm=LogNorm())

        cbar = plt.colorbar(ticks=[10.,100.,1000.])
        cbar.ax.set_yticklabels(['10','100','1000'])

        plt.title(eventfile)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.show() 
    
    
# Function to plot GTI (Good Time Interval) filtered images
def gti_filtered_image(eventfile, out_clean_evtFile, coordinates, expression):
    
    plt.figure(figsize=(20,7))

    pl=1
    
    hdu_list = fits.open(eventfile, memmap=True)
    evt_data = Table(hdu_list[1].data)
    prihdu   = hdu_list[1].header
    print("Events in event file" + " " + eventfile + ": " + str(len(evt_data)) + "\n")
    
    gti_hdu_list = fits.open(out_clean_evtFile, memmap=True)
    gti_evt_data = Table(gti_hdu_list[1].data)
    print("Events in GTI clean event file" + " " + out_clean_evtFile + ": " + str(len(gti_evt_data)) + "\n")
    
    if coordinates == 'DET':
        x_data = evt_data['DETX']
        y_data = evt_data['DETY']
        x_label = 'DETX'
        y_label = 'DETY'
    else:
        x_data = evt_data['X']
        y_data = evt_data['Y']
        x_label = 'X'
        y_label = 'Y'


    # Create Events image
    
    xmax=np.amax(x_data)
    xmin=np.amin(x_data)
    xmid=abs(xmax-xmin)/2.+xmin
    ymax=np.amax(y_data)
    ymin=np.amin(y_data)
    xbin_size=80
    ybin_size=80
    NBINS = (int(abs(xmax-xmin)/xbin_size),int(abs(ymax-ymin)/ybin_size))
    
    plt.subplot(1, 2, pl)
    
    img_zero_mpl = plt.hist2d(x_data, y_data, NBINS, cmap='GnBu', norm=LogNorm())
    
    cbar = plt.colorbar(ticks=[10.,100.,1000.])
    cbar.ax.set_yticklabels(['10','100','1000'])
    
    plt.title(eventfile)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    pl=pl+1
    
    # Create Clean Events image
    if coordinates == 'DET':
        gti_x_data = gti_evt_data['DETX']
        gti_y_data = gti_evt_data['DETY']
        gti_x_label = 'DETX'
        gti_y_label = 'DETY'
    else:
        gti_x_data = gti_evt_data['X']
        gti_y_data = gti_evt_data['Y']
        gti_x_label = 'X'
        gti_y_label = 'Y'


    xmax=np.amax(gti_x_data)
    xmin=np.amin(gti_x_data)
    xmid=abs(xmax-xmin)/2.+xmin
    ymax=np.amax(gti_y_data)
    ymin=np.amin(gti_y_data)
    xbin_size=80
    ybin_size=80
    NBINS = (int(abs(xmax-xmin)/xbin_size),int(abs(ymax-ymin)/ybin_size))
    
    plt.subplot(1, 2, pl)
    
    img_zero_mpl = plt.hist2d(gti_x_data, gti_y_data, NBINS, cmap='GnBu', norm=LogNorm())
    
    cbar = plt.colorbar(ticks=[10.,100.,1000.])
    cbar.ax.set_yticklabels(['10','100','1000'])
    
    plt.title(out_clean_evtFile)
    plt.xlabel(gti_x_label)
    plt.ylabel(gti_y_label)
    
    plt.text(xmid, ymin-0.1*abs(ymax-ymin), expression, ha='center')
    
    pl=pl+1
    
    gti_hdu_list.close()
    hdu_list.close()   
    