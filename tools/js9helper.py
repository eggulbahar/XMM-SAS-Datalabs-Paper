"""
Author: Esin G. Gulbahar
Last Update: 13/12/2024
This code provides a tool for visualizing and analyzing astronomical images using the jpyjs9 and astropy libraries. It includes two main functions:

visualise: This function displays a FITS image in the JS9 interface, allowing users to adjust the scale, colormap, contrast, and bias for better visualization.

getRegions: This function retrieves and processes region data from the JS9 interface. It identifies and distinguishes between "source" and "background" regions, extracting coordinates and region details (like radius, RA/Dec) for various coordinate systems (e.g., FK5, physical, image, ecliptic, galactic). The results are returned in a dictionary, with values such as the coordinates and sizes of the identified regions. It saves the region coordinate file.
"""

import jpyjs9
from astropy.io import fits
import numpy as np
import json
import os

def visualise(my_js9, out_IMFile, scale= 'log', colormap = 'bb', contrast = 2, bias = 0.4):
    hdul = fits.open(out_IMFile)
    my_js9.SetFITS(hdul)
    my_js9.SetColormap(colormap, contrast, bias)
    my_js9.SetScale(scale)
    
def getRegions(my_js9, filename="regions_data.json"):
    my_regions = my_js9.GetRegions()
    
    source_reg = None  
    bkg_reg = None
    
    print(my_regions)
    source_found = False
    
    result = {'x_source': None, 'y_source': None, 'r_source': None, 'r2_source': None, 'ra_source':None, 'dec_source':None, 'radius_src':None, 'radius1_src':None, 'radius2_src': None, 'ra_bkg': None, 'dec_bkg': None, 'radius_bkg':None, 'x_bkg': None, 'y_bkg': None, 'r_bkg': None, 'r2_bkg': None}
    for region in my_regions:
        
        # Check if any tag in the 'tags' list contains 'source' or 'src' (case-insensitive)
        if any(tag.lower() in ['source', 'src'] for tag in region['tags']):
            source_found = True
            print("Found source or src tag")
            source = region
            source_reg = source['imstr']
            source_reg2 = None
            ra = source['ra']
            dec = source['dec']
            
            
            coord_sys = None
            if 'wcsconfig' in region and 'wcssys' in region['wcsconfig']:
                if region['wcsconfig']['wcssys'] in ['physical', 'image', 'ecliptic', 'galactic', 'FK5']:
                    coord_sys = region['wcsconfig']['wcssys']
            if coord_sys:
                print(f"Determined the coordinate system: {coord_sys}")
                
            if 'wcsconfig' in region and 'wcsstr' in region['wcsconfig']:
                source_reg2 = region['wcsconfig']['wcsstr']
                
            if coord_sys == 'FK5':
                if region['shape'] == 'circle':
                    result['ra_source'] = ra
                    result['dec_source'] = dec
                    result['radius_src'] = str(float(source_reg2.split(",")[2].replace('")',''))/3600)
                    print("The coordinates of the selected source region are: \n")
                    print("  ra_source = ", result['ra_source'], "(FK5)")
                    print("  dec_source = ", result['dec_source'], "(FK5)")
                    print("  radius_src = ", result['radius_src'], "(FK5, converted from arcsec) \n")
                    result['x_source'] = source_reg.split(",")[0].replace('circle(','')
                    result['y_source'] = source_reg.split(",")[1].replace('','')
                    result['r_source'] = source_reg.split(",")[2].replace(')','')
                    print("  x_source = ", result['x_source'], "(physical)")
                    print("  y_source = ", result['y_source'], "(physical)")
                    print("  r_source = ", result['r_source'], "(physical) \n")
                elif region['shape'] == 'annulus':
                    result['ra_source'] = ra
                    result['dec_source'] = dec
                    result['radius1_src'] = str(float(source_reg2.split(",")[2].replace('"',''))/3600)
                    result['radius2_src'] = str(float(source_reg2.split(",")[3].replace('")',''))/3600)
                    print("The coordinates of the selected source region are: \n")
                    print("  ra_source = ", result['ra_source'], "(FK5)")
                    print("  dec_source = ", result['dec_source'], "(FK5)")
                    print("  radius_src = ", result['radius_src'], "(FK5, converted from arcsec) \n")
                    print("  radius2_src = ", result['radius2_src'], "(FK5, converted from arcsec) \n")
                    result['x_source'] = source_reg.split(",")[0].replace('annulus(','')
                    result['y_source'] = source_reg.split(",")[1].replace('','')
                    result['r_source'] = source_reg.split(",")[2].replace('"','')
                    result['r2_source'] = source_reg.split(",")[3].replace('")','')
                    print("  x_source = ", result['x_source'], "(physical)")
                    print("  y_source = ", result['y_source'], "(physical)")
                    print("  r_source = ", result['r_source'], "(physical)")
                    print("  r2_source = ", result['r2_source'], "(physical)\n")
            
            elif coord_sys == 'physical':
                
                if region['shape'] == 'circle':
                    
                    result['x_source'] = source_reg.split(",")[0].replace('circle(','')
                    result['y_source'] = source_reg.split(",")[1].replace('','')
                    result['r_source'] = source_reg.split(",")[2].replace(')','')
                    print("The coordinates of the selected source region are: \n")
                    print("  x_source = ", result['x_source'], "(physical)")
                    print("  y_source = ", result['y_source'], "(physical)")
                    print("  r_source = ", result['r_source'], "(physical) \n")
                elif region['shape'] == 'annulus':
                    result['x_source'] = source_reg.split(",")[0].replace('annulus(','')
                    result['y_source'] = source_reg.split(",")[1].replace('','')
                    result['r_source'] = source_reg.split(",")[2].replace('"','')
                    result['r2_source'] = source_reg.split(",")[3].replace('")','')
                    print("The coordinates of the selected region are: \n")
                    print("  x_source = ", result['x_source'], "(physical)")
                    print("  y_source = ", result['y_source'], "(physical)")
                    print("  r_source = ", result['r_source'], "(physical)")
                    print("  r2_source = ", result['r2_source'], "(physical)\n")
            elif coord_sys == 'image':
                if region['shape'] == 'circle':
                    result['x_source'] = source_reg.split(",")[0].replace('circle(','')
                    result['y_source'] = source_reg.split(",")[1].replace('','')
                    result['r_source'] = source_reg.split(",")[2].replace(')','')
                    print("The coordinates of the selected source region are: \n")
                    print("  x_source = ", result['x_source'], "(image)")
                    print("  y_source = ", result['y_source'], "(image)")
                    print("  r_source = ", result['r_source'], "(image) \n")
                elif region['shape'] == 'annulus':
                    result['x_source'] = source_reg.split(",")[0].replace('annulus(','')
                    result['y_source'] = source_reg.split(",")[1].replace('','')
                    result['r_source'] = source_reg.split(",")[2].replace('"','')
                    result['r2_source'] = source_reg.split(",")[3].replace('")','')
                    print("The coordinates of the selected region are: \n")
                    print("  x_source = ", result['x_source'], "(image)")
                    print("  y_source = ", result['y_source'], "(image)")
                    print("  r_source = ", result['r_source'], "(image)")
                    print("  r2_source = ", result['r2_source'], "(image)\n")
            elif coord_sys ==  'ecliptic':
                if region['shape'] == 'circle':
                    result['x_source'] = source_reg2.split(",")[0].replace('circle(','')
                    result['y_source'] = source_reg2.split(",")[1].replace('','')
                    result['r_source'] = source_reg2.split(",")[2].replace(')','')
                    print("The coordinates of the selected source region are: \n")
                    print("  x_source = ", result['x_source'], "(ecliptic)")
                    print("  y_source = ", result['y_source'], "(ecliptic)")
                    print("  r_source = ", result['r_source'], "(ecliptic) \n")
                elif region['shape'] == 'annulus':
                    result['x_source'] = source_reg2.split(",")[0].replace('annulus(','')
                    result['y_source'] = source_reg2.split(",")[1].replace('','')
                    result['r_source'] = source_reg2.split(",")[2].replace('"','')
                    result['r2_source'] = source_reg2.split(",")[3].replace('")','')
                    print("The coordinates of the selected region are: \n")
                    print("  x_source = ", result['x_source'], "(ecliptic)")
                    print("  y_source = ", result['y_source'], "(ecliptic)")
                    print("  r_source = ", result['r_source'], "(ecliptic)")
                    print("  r2_source = ", result['r2_source'], "(ecliptic)\n")
            elif coord_sys == 'galactic':
                if region['shape'] == 'circle':
                    result['x_source'] = source_reg2.split(",")[0].replace('circle(','')
                    result['y_source'] = source_reg2.split(",")[1].replace('','')
                    result['r_source'] = source_reg2.split(",")[2].replace(')','')
                    print("The coordinates of the selected source region are: \n")
                    print("  x_source = ", result['x_source'], "(galactic)")
                    print("  y_source = ", result['y_source'], "(galactic)")
                    print("  r_source = ", result['r_source'], "(galactic) \n")
                elif region['shape'] == 'annulus':
                    result['x_source'] = source_reg2.split(",")[0].replace('annulus(','')
                    result['y_source'] = source_reg2.split(",")[1].replace('','')
                    result['r_source'] = source_reg2.split(",")[2].replace('"','')
                    result['r2_source'] = source_reg2.split(",")[3].replace('")','')
                    print("The coordinates of the selected region are: \n")
                    print("  x_source = ", result['x_source'], "(galactic)")
                    print("  y_source = ", result['y_source'], "(galactic)")
                    print("  r_source = ", result['r_source'], "(galactic)")
                    print("  r2_source = ", result['r2_source'], "(galactic)\n")
    

           
            
        if any(tag.lower() in ['back', 'background', 'bkg', 'bg'] for tag in region['tags']):
            source_found = True
            print("Found back, background, bg or bkg tag")
            background = region
            bkg_reg = background['imstr']
            bkg_reg2 = None
            ra = background['ra']
            dec = background['dec']
            
            coord_sys = None
            if 'wcsconfig' in region and 'wcssys' in region['wcsconfig']:
                if region['wcsconfig']['wcssys'] in ['physical', 'image', 'ecliptic', 'galactic', 'FK5']:
                    coord_sys = region['wcsconfig']['wcssys']
            if coord_sys:
                print(f"Determined the coordinate system: {coord_sys}")
                
            if 'wcsconfig' in region and 'wcsstr' in region['wcsconfig']:
                bkg_reg2 = region['wcsconfig']['wcsstr']
            
            if coord_sys == 'FK5':
                if region['shape'] == 'circle':
                    result['ra_bkg'] = ra
                    result['dec_bkg'] = dec
                    result['radius_bkg'] = str(float(bkg_reg2.split(",")[2].replace('")',''))/3600)
                    print("The coordinates of the selected source region are: \n")
                    print("  ra_bkg = ", result['ra_bkg'], "(FK5)")
                    print("  dec_bkg = ", result['dec_bkg'], "(FK5)")
                    print("  radius_bkg = ", result['radius_bkg'], "(FK5, converted from arcsec) \n")
                    result['x_bkg'] = bkg_reg.split(",")[0].replace('circle(','')
                    result['y_bkg'] = bkg_reg.split(",")[1].replace('','')
                    result['r_bkg'] = bkg_reg.split(",")[2].replace(')','')
                    print("  x_bkg = ", result['x_bkg'], "(physical)")
                    print("  y_bkg = ", result['y_bkg'], "(physical)")
                    print("  r_bkg = ", result['r_bkg'], "(physical)\n")
                elif region['shape'] == 'annulus':
                    result['ra_bkg'] = ra
                    result['dec_bkg'] = dec
                    result['radius1_bkg'] = str(float(bkg_reg2.split(",")[2].replace('"',''))/3600)
                    result['radius2_bkg'] = str(float(bkg_reg2.split(",")[3].replace('")',''))/3600)
                    print("The coordinates of the selected source region are: \n")
                    print("  ra_bkg = ", result['ra_bkg'], "(FK5)")
                    print("  dec_bkg = ", result['dec_bkg'], "(FK5)")
                    print("  radius1_bkg = ", result['radius1_bkg'], "(FK5, converted from arcsec) \n")
                    print("  radius2_bkg = ", result['radius2_bkg'], "(FK5, converted from arcsec) \n")
                    result['x_bkg'] = bkg_reg.split(",")[0].replace('annulus(','')
                    result['y_bkg'] = bkg_reg.split(",")[1].replace('','')
                    result['r_bkg'] = bkg_reg.split(",")[2].replace('','')
                    result['r2_bkg'] = bkg_reg.split(",")[3].replace(')','')
                    print("  x_bkg = ", result['x_bkg'], "(physical)")
                    print("  y_bkg = ", result['y_bkg'], "(physical)")
                    print("  r_bkg = ", result['r_bkg'], "(physical)")
                    print("  r2_bkg = ", result['r2_bkg'], "(physical)\n")
            
            elif coord_sys == 'physical':
                if region['shape'] == 'circle':
                    result['x_bkg'] = bkg_reg.split(",")[0].replace('circle(','')
                    result['y_bkg'] = bkg_reg.split(",")[1].replace('','')
                    result['r_bkg'] = bkg_reg.split(",")[2].replace(')','')
                    print("The coordinates of the selected background region are: \n")
                    print("  x_bkg = ", result['x_bkg'], "(physical)")
                    print("  y_bkg = ", result['y_bkg'], "(physical)")
                    print("  r_bkg = ", result['r_bkg'], "(physical)\n")
                elif region['shape'] == 'annulus':
                    result['x_bkg'] = bkg_reg.split(",")[0].replace('annulus(','')
                    result['y_bkg'] = bkg_reg.split(",")[1].replace('','')
                    result['r_bkg'] = bkg_reg.split(",")[2].replace('','')
                    result['r2_bkg'] = bkg_reg.split(",")[3].replace(')','')
                    print("The coordinates of the selected region are: \n")
                    print("  x_bkg = ", result['x_bkg'], "(physical)")
                    print("  y_bkg = ", result['y_bkg'], "(physical)")
                    print("  r_bkg = ", result['r_bkg'], "(physical)")
                    print("  r2_bkg = ", result['r2_bkg'], "(physical)\n")
            elif coord_sys == 'image':
                if region['shape'] == 'circle':
                    result['x_bkg'] = bkg_reg.split(",")[0].replace('circle(','')
                    result['y_bkg'] = bkg_reg.split(",")[1].replace('','')
                    result['r_bkg'] = bkg_reg.split(",")[2].replace(')','')
                    print("The coordinates of the selected background region are: \n")
                    print("  x_bkg = ", result['x_bkg'], "(image)")
                    print("  y_bkg = ", result['y_bkg'], "(image)")
                    print("  r_bkg = ", result['r_bkg'], "(image)\n")
                elif region['shape'] == 'annulus':
                    result['x_bkg'] = bkg_reg.split(",")[0].replace('annulus(','')
                    result['y_bkg'] = bkg_reg.split(",")[1].replace('','')
                    result['r_bkg'] = bkg_reg.split(",")[2].replace('','')
                    result['r2_bkg'] = bkg_reg.split(",")[3].replace(')','')
                    print("The coordinates of the selected region are: \n")
                    print("  x_bkg = ", result['x_bkg'], "(image)")
                    print("  y_bkg = ", result['y_bkg'], "(image)")
                    print("  r_bkg = ", result['r_bkg'], "(image)")
                    print("  r2_bkg = ", result['r2_bkg'], "(image)\n")
            elif coord_sys == 'ecliptic':
                if region['shape'] == 'circle':
                    result['x_bkg'] = bkg_reg2.split(",")[0].replace('circle(','')
                    result['y_bkg'] = bkg_reg2.split(",")[1].replace('','')
                    result['r_bkg'] = bkg_reg2.split(",")[2].replace(')','')
                    print("The coordinates of the selected background region are: \n")
                    print("  x_bkg = ", result['x_bkg'], "(ecliptic)")
                    print("  y_bkg = ", result['y_bkg'], "(ecliptic)")
                    print("  r_bkg = ", result['r_bkg'], "(ecliptic)\n")
                elif region['shape'] == 'annulus':
                    result['x_bkg'] = bkg_reg2.split(",")[0].replace('annulus(','')
                    result['y_bkg'] = bkg_reg2.split(",")[1].replace('','')
                    result['r_bkg'] = bkg_reg2.split(",")[2].replace('','')
                    result['r2_bkg'] = bkg_reg2.split(",")[3].replace(')','')
                    print("The coordinates of the selected region are: \n")
                    print("  x_bkg = ", result['x_bkg'], "(ecliptic)")
                    print("  y_bkg = ", result['y_bkg'], "(ecliptic)")
                    print("  r_bkg = ", result['r_bkg'], "(ecliptic)")
                    print("  r2_bkg = ", result['r2_bkg'], "(ecliptic)\n")
            elif coord_sys == 'galactic':
                if region['shape'] == 'circle':
                    result['x_bkg'] = bkg_reg2.split(",")[0].replace('circle(','')
                    result['y_bkg'] = bkg_reg2.split(",")[1].replace('','')
                    result['r_bkg'] = bkg_reg2.split(",")[2].replace(')','')
                    print("The coordinates of the selected background region are: \n")
                    print("  x_bkg = ", result['x_bkg'], "(galactic)")
                    print("  y_bkg = ", result['y_bkg'], "(galactic)")
                    print("  r_bkg = ", result['r_bkg'], "(galactic)\n")
                elif region['shape'] == 'annulus':
                    result['x_bkg'] = bkg_reg2.split(",")[0].replace('annulus(','')
                    result['y_bkg'] = bkg_reg2.split(",")[1].replace('','')
                    result['r_bkg'] = bkg_reg2.split(",")[2].replace('','')
                    result['r2_bkg'] = bkg_reg2.split(",")[3].replace(')','')
                    print("The coordinates of the selected region are: \n")
                    print("  x_bkg = ", result['x_bkg'], "(galactic)")
                    print("  y_bkg = ", result['y_bkg'], "(galactic)")
                    print("  r_bkg = ", result['r_bkg'], "(galactic)")
                    print("  r2_bkg = ", result['r2_bkg'], "(galactic)\n")
    # Save results to a file
    if os.path.exists(filename):
        with open(filename, "r") as file:
            try:
                existing_data = json.load(file)
            except json.JSONDecodeError:
                existing_data = []
    else:
        existing_data = []
    
    existing_data.append(result)  # Append new result to existing data
    
    with open(filename, "w") as file:
        json.dump(existing_data, file, indent=4)
    
    print(f"Regions data saved to {filename}")
    return result

    
    
    if not source_found:
        print("No source or background region selected.")
            
        
    
    

   
