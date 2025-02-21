"""
Author: Esin G. Gulbahar
Last updated: 13/12/2024

This Python script provides tools for analyzing and visualizing X-ray astronomy light curves, with a focus on data from FITS files. It integrates multiple libraries to create flexible and interactive visualizations and ensures compatibility with XMM-Newton datasets. 

Features include:
1. plotLCmatplotlib:
   - Visualizes light curves using Matplotlib.
   - Supports error bars and optional threshold lines.
   - Allows customization of axis labels, limits, and titles.

2. plotLC:
   - Uses Plotly for interactive light curve visualizations.
   - Supports threshold display, dynamic legend updates, and error bar plotting.
   - Offers an intuitive interface for exploring light curve data.

3. read_lightcurve:
   - Converts XMM-Newton light curve FITS files into generic `LightCurve` objects.
   - Handles data preprocessing, such as unit corrections, metadata updates, and time normalization.

4. lcviz:
   - Facilitates advanced visualization of light curves using LCviz.
   - Supports multiple light curves with customizable labels.

Additional Features:
- Interactive plotting using Plotly.
- Error handling for missing data and units.
- Metadata extraction from FITS headers to enrich analysis.
- Ensures compatibility with XMM-Newton data conventions.
"""
import os.path
from os import path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = 'notebook'
import plotly.offline as pyo
pyo.init_notebook_mode(connected=True)
from lcviz import *
from lightkurve import * 
import lightkurve 
import logging
import warnings
from copy import deepcopy
from astropy.time import Time
from astropy.units import UnitsWarning


def plotVelaX1LC(fileNames, names, threshold=None, figname="lightcurve.png", yLog=False, connect_points=False):
    seconds_in_day = 86400  # Number of seconds in a day

    # Create the main plot and axis
    fig, ax1 = plt.subplots(figsize=(12, 6))
    ax1.xaxis.get_offset_text().set_visible(False)

    for fileName, name in zip(fileNames, names):
        if fileName != "NOT FOUND":
            fitsFile = fits.open(fileName)
            prihdu = fitsFile[1].header

            # Extract threshold value from header if available
            if 'CUTVAL' in prihdu:
                threshold = prihdu['CUTVAL']

            cols = fitsFile[1].columns
            colName = None

            # Determine the correct column name for RATE or COUNTS
            for i, x in enumerate(cols.names):
                if "RATE" in x:
                    colName = cols.names[i]
                if "COUNTS" in x:
                    colName = cols.names[i]

            data = fitsFile[1].data
            xdata = data.field('TIME')  # Extract the time column
            ydata = data.field(colName)
            
            mjd_reference = 50814  # Reference MJD
            mjd_target = 58607  # Target MJD
            seconds_in_day = 86400
            # Calculate the reference time difference in seconds
            time_difference = (mjd_target - mjd_reference) * seconds_in_day
            # Convert TT times to relative MJD
            relative_seconds = xdata - time_difference
            relative_mjd = mjd_target + (relative_seconds / seconds_in_day)

            # Convert time from seconds to days since MJD reference
            xdata_days = relative_mjd-58607
            
            

            # Calculate orbital phase with respect to T90
            orbital_phase = ((relative_mjd - 52974.001) / 8.964357) % 1

            xmax = np.amax(xdata_days)
            xmin = np.amin(xdata_days)

            if connect_points:
                ax1.plot(xdata_days, ydata, label=name, linestyle='-', marker='.')
            else:
                ax1.plot(xdata_days, ydata, label=name, linestyle="", marker=".")

            # Set plot labels and titles
            if colName == 'RATE':
                ax1.set_title("XMM EPIC-pn (0.5-10 keV)")
                ax1.set_xlabel(f"Time (days since MJD 58607)")
                ax1.set_ylabel("Cts/s")
            else:
                ax1.set_title("XMM EPIC-pn (0.5-10 keV)")
                ax1.set_xlabel(f"Time (days since MJD 58607)")
                ax1.set_ylabel("Counts")

            # Add a threshold line if specified
            if threshold is not None and threshold != 'None':
                if colName == 'COUNTS':
                    threshold = float(threshold) * 100.

                y2data = [threshold] * len(xdata_days)
                ax1.plot(xdata_days, y2data, linestyle='--', color='red')
                ax1.text(xmin + 0.1 * (xmax - xmin), threshold + 0.01 * threshold,
                         str(threshold) + " cts/sec", ha='center', color='red')

            fitsFile.close()
        else:
            print("File not found " + fileName + "\n")
    
    if yLog:
        ax1.set_yscale('log')


    ax2 = ax1.secondary_xaxis('top')
    ax2.set_xlabel("Orbital Phase")

    # Set the secondary x-axis tick labels as orbital phases
    ax2.set_xticks(ax1.get_xticks())  # Use the same x-ticks as the main axis
    ax2.set_xticklabels([f"{((x +58607 - 52974.001) / 8.964357) % 1:.2f}" for x in ax1.get_xticks()])
    
    
    plt.legend()
    plt.savefig(figname)
    plt.show()



def plotLC(fileNames, names, threshold=None):
    """
    Plots light curves using Plotly from a list of FITS files. Supports optional binning and threshold display.

    Parameters:
    - fileNames: list of str
        List of FITS file paths containing the light curve data.
    - names: list of str
        Corresponding names for each light curve to be used in the legend.
    - threshold: float, optional
        Threshold value to plot as a horizontal line (default: None).

    Returns:
    - None
    """

    # Create an empty Plotly figure
    fig = go.Figure()
    
    # Iterate over each FITS file and its corresponding name
    for fileName, name in zip(fileNames, names):
        if fileName != "NOT FOUND":  # Ensure the file exists
            fitsFile = fits.open(fileName)  # Open the FITS file
            prihdu = fitsFile[1].header  # Access the primary header of the first extension

            # Check if a threshold value is stored in the header
            if 'CUTVAL' in prihdu:
                threshold = prihdu['CUTVAL']

            cols = fitsFile[1].columns  # Access the column information of the data table
            colName = None
            errName = None

            # Identify the appropriate column for light curve data (RATE or COUNTS)
            for i, x in enumerate(cols.names):
                if "RATE" in x:
                    colName = cols.names[i]
                if "COUNTS" in x:
                    colName = cols.names[i]
                if "ERR" in x or "ERROR" in x:
                    errName = cols.names[i]  # Identify error column

            data = fitsFile[1].data  # Extract the data table

            # Extract time and light curve values, normalizing time to start from zero
            xdata = data.field('TIME') - min(data.field('TIME'))
            ydata = data.field(colName)

            # Extract y-errors if available
            yerr = data.field(errName) if errName else None

            # Determine the range of the x-axis
            xmax = np.amax(xdata)
            xmin = np.amin(xdata)

            # Set y-axis label based on the column type
            ylabel = "Cts/s" if colName == 'RATE' else "Counts"

            # Plot threshold line if provided
            if threshold is not None and threshold != 'None':
                if colName == 'COUNTS':  # Convert threshold to equivalent counts if necessary
                    threshold = float(threshold) * 100.

                y2data = [threshold] * len(xdata)  # Create a constant threshold line
                fig.add_trace(go.Scatter(x=xdata, y=y2data, name="Threshold (cts/sec)", mode='lines'))

            # Plot light curve with error bars if yerr is available
            if yerr is not None:
                fig.add_trace(go.Scatter(x=xdata, y=ydata, mode='lines', name="Lightcurve of " + name, 
                                         error_y=dict(type='data', array=yerr, visible=True)))
            else:
                fig.add_trace(go.Scatter(x=xdata, y=ydata, mode='lines', name="Lightcurve of " + name))

            fitsFile.close()  # Close the FITS file

        else:
            # Print a warning if the file is not found
            print(f"File not found: {fileName}\n")

    # Update the overall plot layout
    fig.update_layout(title='Lightcurve',
                      xaxis_title='Time (s)',
                      yaxis_title=ylabel)
    
    # Display the plot
    fig.show()

 
def read_lightcurve(
    fileName,
    time_column="time",
    flux_column="rate",
    flux_err_column="error",
    quality_column="quality",
    cadenceno_column="cadenceno",
    centroid_col_column="mom_centr1",
    centroid_row_column="mom_centr2",
    time_format="MET",  # Default time format as MET
    ext=1,
):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        """Generic helper function to convert XMM-Newton light curve file
        into a generic `LightCurve` object.
        """

        # Open the FITS file
        if isinstance(fileName, fits.HDUList):
            hdulist = fileName  # Allow HDUList to be passed
        else:
            with fits.open(fileName) as hdulist:
                hdulist = deepcopy(hdulist)

        # Check if the requested extension exists
        if isinstance(ext, str):
            validate_method(ext, supported_methods=[hdu.name.lower() for hdu in hdulist])

        # Read the data table
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UnitsWarning)
            tab = Table.read(fileName, format="fits")

        # Check and update metadata
        tab.meta.update(hdulist[1].header)
        tab.meta = {k: v for k, v in tab.meta.items()}

        # Ensure columns have correct units
        for colname in tab.colnames:
            unitstr = str(tab[colname].unit)
            if unitstr == "e-/s":
                tab[colname].unit = "electron/s"
            elif unitstr == "ct/s":
                tab[colname].unit = "cts/s"
            elif unitstr.lower() == "unitless":
                tab[colname].unit = ""

            tab.rename_column(colname, colname.lower())

        # Rename TIME column to lowercase "time" if necessary
        if time_column == "TIME" and "time" not in tab.columns and "TIME" in tab.colnames:
            tab.rename_column("TIME", "time")

        # Remove rows with NaN time values
        nans = np.isnan(tab["time"].data)
        if np.any(nans):
            tab = tab[~nans]

        # Verify the time units are in seconds
        time_unit = hdulist[ext].header.get("TIMEUNIT", "s").lower()
        if time_unit != "s":
            raise ValueError(f"Unexpected time unit '{time_unit}'. Expected 's'.")

        # Prepare the time column as MET (Unix format)
        if time_format == "MET":
            # Use 'unix' format, suitable for Mission Elapsed Time (MET)
            time = Time(tab["time"].data, scale=hdulist[ext].header.get("TIMESYS", "tdb").lower(), format="unix")
        else:
            time = Time(
                tab["time"].data,
                scale=hdulist[ext].header.get("TIMESYS", "tdb").lower(),
                format="unix",  # If another time format is desired
            )

        # For backwards compatibility, ensure standard columns exist
        if flux_column not in tab.columns:
            alt_flux_column = f"{flux_column}_alt"
            if alt_flux_column in tab.columns:
                tab.add_column(tab[alt_flux_column], name=flux_column)
            else:
                raise KeyError(f"Column '{flux_column}' not found in the table.")

        if flux_err_column not in tab.columns:
            raise KeyError(f"Column '{flux_err_column}' not found in the table.")

        if "quality" not in tab.columns and quality_column in tab.columns:
            tab.add_column(tab[quality_column], name="quality", index=2)

        if "cadenceno" not in tab.columns and cadenceno_column in tab.columns:
            tab.add_column(tab[cadenceno_column], name="cadenceno", index=3)

        if "centroid_col" not in tab.columns and centroid_col_column in tab.columns:
            tab.add_column(tab[centroid_col_column], name="centroid_col", index=4)

        if "centroid_row" not in tab.columns and centroid_row_column in tab.columns:
            tab.add_column(tab[centroid_row_column], name="centroid_row", index=5)

        # Update metadata
        tab.meta["LABEL"] = hdulist[0].header.get("OBJECT")
        tab.meta["MISSION"] = hdulist[0].header.get(
            "MISSION", hdulist[0].header.get("TELESCOP")
        )
        tab.meta["RA"] = hdulist[0].header.get("RA_OBJ")
        tab.meta["DEC"] = hdulist[0].header.get("DEC_OBJ")
        tab.meta["FILENAME"] = fileName
        tab.meta["FLUX_ORIGIN"] = flux_column

        # Return the LightCurve object, but without passing time explicitly since it's in the table
        return LightCurve(data=tab, flux=tab[flux_column], flux_err=tab[flux_err_column])




def lcviz(xmm_lightcurves, labels=None, threshold = None):
    """
    Visualize multiple light curves using LCviz.

    Parameters:
        xmm_lightcurves (list): List of light curve files converted to an LightCurve object.
        labels (list): Optional list of labels for each light curve. 
                       If None, generic labels will be generated.
    """
    # Create LCviz viewer
    lcviz = LCviz()
    
    # Generate default labels if none are provided
    if labels is None:
        labels = [f"Light Curve {i+1}" for i in range(len(xmm_lightcurves))]
    
    # Ensure the number of labels matches the number of light curves
    if len(labels) != len(xmm_lightcurves):
        raise ValueError("The number of labels must match the number of light curves.")
    
    # Load and add each light curve to LCviz
    for i, xmm_lightcurve in enumerate(xmm_lightcurves):
        lc_src = LightCurve(xmm_lightcurve)
        lcviz.load_data(lc_src, data_label=labels[i])
        
    if threshold is not None:
        # Create a synthetic LightCurve with constant flux
        time = xmm_lightcurves[0].time  # Use the same time grid as the first light curve
        flux = np.full_like(time, threshold)  # Constant flux values at y_line
        flux_err = np.zeros_like(time)  # No errors for the horizontal line
        
        # Create a synthetic light curve
        constant_lc = LightCurve(data={'time': time, 'flux': flux, 'flux_err': flux_err})
        
        # Add it to the viewer
        lcviz.load_data(constant_lc, data_label=f"Threshold Line at {threshold}")
        plotopt=lcviz.plugins['Plot Options']
        plotopt.line_width = 3
        plotopt.line_visible = True
    
    # Customize axis labels
    #lcviz.viewers['flux-vs-time']._obj.state.x_axislabel = 'Seconds'
    lcviz.viewers['flux-vs-time']._obj.state.y_axislabel = 'Count rate'
    po = lcviz.plugins['Plot Options']
    #print(f"viewer choices: {po.viewer.choices}")
    #po.viewer = po.viewer.choices[0]
    #print(f"layer choices: {po.layer.choices}")
    #po.layer = po.layer.choices[0]
    #po.line_color = 'red'
    po.line_width = 2
    po.marker_size = 4
    po.line_visible = True
    # Show the visualization
    lcviz.show()

"""
def inspect_time_format(file_name):
    with fits.open(file_name) as hdul:
        # Access the primary header
        primary_header = hdul[0].header
        # Access the extension header containing the light curve data
        lc_header = hdul[1].header

        # Print relevant time-related keywords
        print("Primary Header:")
        for keyword in ["TIMESYS", "MJDREF", "BJDREF", "TIMEUNIT", "TIMEZERO"]:
            if keyword in primary_header:
                print(f"{keyword}: {primary_header[keyword]}")

        print("\nLight Curve Header:")
        for keyword in ["TIMESYS", "MJDREF", "BJDREF", "TIMEUNIT", "TIMEZERO", "TSTART", "TSTOP", "TIMEDEL"]:
            if keyword in lc_header:
                print(f"{keyword}: {lc_header[keyword]}")
"""