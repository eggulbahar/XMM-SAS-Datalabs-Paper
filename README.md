# XMM-Newton SAS in ESA Datalabs
### *A New Browser Based Interactive Environment for XMM-Newton Analysis*

This is a repository with code and data analysis notebook for the X-ray binary Vela X-1 presented in the paper: .

Check [ReadMeTools.md](ReadMeTools.md) for details on the Python functions.
Check [ReadMeScripts.md](ReadMeScripts.md) for details on the Python scripts.

## Folder structure of the repository (XMM-SAS-Datalabs-Paper):

* [CaseStudy.ipynb](CaseStudy.ipynb): Jupyter Lab Notebook with the pySAS data extraction and visualisation on Vela X-1.
* [tools](tools): Python utility functions needed for plotting and visualisation.
* [scripts](scripts): Python scripts to run looped SAS tasks to extract data.
* [data](data): folder with all the output data.

## Pre-requisites
If running the Notebook on [ESA Datalabs](https://datalabs.esa.int/) inside the XMM-SAS datalab, no pre-requisites are required.

To run locally you need python 3.10 or above and standard packages like `astropy`, `numpy`, `scipy`, `plotly`, `lcviz` and [HEASOFT](https://heasarc.gsfc.nasa.gov/lheasoft/) built from source in order to have `pyXspec` available. You will need SAS version 21.

---

*Author: Esin G. Gulbahar*

*Last Updated: 24/02/2025*