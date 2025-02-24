# Python Tools for Visualisation

This directory, named **"tools"**, is a collection of Python functions designed to simplify data visualisation. These files include interactive plotting and image visualization functions tailored for XMM-Newton data. The functions can be customised by copying them to a personal directory.

---

## **Contents**

### **1. [js9helper.py](js9helper.py)**  
A utility script to simplify the use of the JS9 tool in a JupyterLab environment.  
- **Functions:**
  - `visualise`: Displays FITS images within the JS9 interface.  
  - `getRegions`: Retrieves and processes region data for further analysis.
 
### **2. [xspecplot.py](xspecplot.py)**  
Generates stacked spectral plots using XSPEC and Matplotlib. It retrieves spectral data (energy values, count rates, and errors) from pyXSPEC, applies logarithmic scaling if specified, and overlays reference lines at specific energies.

### **3. [plotLC.py](plotLC.py)**  
This Python script provides tools for analysing and visualising X-ray astronomy light curves from FITS files, integrating Matplotlib, Plotly, Astropy, and LCviz for customizable and interactive plots, preprocessing XMM-Newton data into LightCurve objects, handling metadata, and supporting advanced visualizations for astrophysical research.
- **Functions:**
  - `plotVelaX1LC`: Plots the light curve of Vela X-1 using Matplotlib. It extracts time and count rate data from FITS files, converts time to days since MJD 58607, and overlays orbital phase information.  
  - `plotLC`: Uses Plotly for interactive light curve visualizations.
  - `read_lightcurve`: Converts XMM-Newton light curve FITS files into generic `LightCurve` objects.
  - `lcviz`: Uses LCviz for interactive light curve visualizations.

---

*Author: Esin G. Gulbahar*

*Last updated: 24/02/2025*