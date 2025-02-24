# Scripts for Looped SAS Tasks

This directory, named **"scripts"**, is a collection of Python scripts designed to simplify looped data extration using pySAS. The scripts can be customised by copying them to a personal directory.

---

## **Contents**

### **1. [energy-resolvedLC.py](energy-resolvedLC.py)**  
This script extracts light curves for the each four bands of interest for further analysis of variability in photon count rates. The final energy-resolved light curves are saved in a directory and appended to a list for further analysis. 
 
### **2. [spectrum-extractor.py](spectrum-extractor.py)**  
This script calculates time intervals for spectral extraction from three orbital phase points, sets up the SAS environment, and iteratively extracts source spectra from specific detector regions based on time filtering. The script then generates response (RMF) and ancillary (ARF) files, applies spectral grouping, and outputs the final grouped spectra.

### **3. [gtiloop.py](gtiloop.py)**  
This script initializes the SAS environment and reads the event file to extract the observation start and end times. It then creates Good Time Interval (GTI) files by iterating over the pulse period in 283.44-second intervals using the `tabgtigen` command.

### **4. [loopgtispectra.py](loopgtispectra.py)** 
This script iterates over GTI files, produced running `gtiloop.py`, to extract spectra while avoiding pile-up regions using `evselect`, then applies background scaling (`backscale`), response matrix generation (`rmfgen`, and ancillary response file creation (`arfgen`). Finally, it groups the spectra using `specgroup` and saves the outputs.

---

*Author: Esin G. Gulbahar*

*Last updated: 24/02/2025*