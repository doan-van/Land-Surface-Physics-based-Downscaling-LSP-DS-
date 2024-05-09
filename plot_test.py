#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 15:31:58 2024

@author: doan
"""




import sys, glob, os
import xarray as xr
from matplotlib import gridspec
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches
from matplotlib.colors import from_levels_and_colors
import matplotlib as mpl
from math import floor
import matplotlib
from matplotlib import patheffects
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


import xarray as xr
import numpy as np
import glob, os, sys
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
from matplotlib import gridspec
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches
from matplotlib.colors import from_levels_and_colors
import matplotlib as mpl
from math import floor
import matplotlib
from matplotlib import patheffects
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature



#==============================================================================
def add_grid(ax,ts,fontsize = 10):
        # grid ----
    import matplotlib.ticker as mticker
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top   = False
    gl.ylabels_right = False
    gl.xlines        = False
    gl.ylines        = False
    gl.ylocator = mticker.FixedLocator(np.arange(-180,180,ts))
    gl.xlocator = mticker.FixedLocator(np.arange(-180,180,ts))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylabel_style = {'size': fontsize, 'color': 'gray', 'rotation': 90, 'va':'center'}
    gl.xlabel_style = {'size': fontsize, 'color': 'gray'} #, 'weight': 'bold'}
    #-----------------------   
    return gl
#==============================================================================



if os.name == 'nt':
    matplotlib.rc('font', family='Arial')
else:  # might need tweaking, must support black triangle for N arrow
    matplotlib.rc('font', family='DejaVu Sans')


def utm_from_lon(lon):
    """
    utm_from_lon - UTM zone for a longitude

    Not right for some polar regions (Norway, Svalbard, Antartica)

    :param float lon: longitude
    :return: UTM zone number
    :rtype: int
    """
    return floor( ( lon + 180 ) / 6) + 1

def scale_bar(ax, proj, length, location=(0.5, 0.05), linewidth=3,
              units='km', m_per_unit=1000, north=False):
    """

    http://stackoverflow.com/a/35705477/1072212
    ax is the axes to draw the scalebar on.
    proj is the projection the axes are in
    location is center of the scalebar in axis coordinates ie. 0.5 is the middle of the plot
    length is the length of the scalebar in km.
    linewidth is the thickness of the scalebar.
    units is the name of the unit
    m_per_unit is the number of meters in a unit
    """
    # find lat/lon center to find best UTM zone
    x0, x1, y0, y1 = ax.get_extent(proj.as_geodetic())
    # Projection in metres
    utm = ccrs.UTM(utm_from_lon((x0+x1)/2))
    # Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(utm)
    # Turn the specified scalebar location into coordinates in metres
    sbcx, sbcy = x0 + (x1 - x0) * location[0], y0 + (y1 - y0) * location[1]
    # Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbcx - length * m_per_unit/2, sbcx + length * m_per_unit/2]
    # buffer for scalebar
    buffer = [patheffects.withStroke(linewidth=5, foreground="w")]
    # Plot the scalebar with buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, path_effects=buffer)
    # buffer for text
    buffer = [patheffects.withStroke(linewidth=3, foreground="w")]
    # Plot the scalebar label
    t0 = ax.text(sbcx, sbcy, str(length) + ' ' + units, transform=utm,
        horizontalalignment='center', verticalalignment='bottom',
        path_effects=buffer, zorder=2)
    left = x0+(x1-x0)*0.05
    # Plot the N arrow
    if north:
        t1 = ax.text(left, sbcy, u'\u25B2\nN', transform=utm,
            horizontalalignment='center', verticalalignment='bottom',
            path_effects=buffer, zorder=2)
        
    # Plot the scalebar without buffer, in case covered by text buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, zorder=3)
    
#==============================================================================






geo = xr.open_dataset('geo_em.d03.nc')
lon, lat = geo.XLONG_M[0].values,geo.XLAT_M[0].values

ifiles = sorted(glob.glob('dataout/extract/TEMP/*.nc'))
for f in ifiles[:-1][-1:]:
    d = xr.open_dataset(f).TEMP
    #d.plot()
                
    temp = d[0] 
    
    
fig   =  plt.figure(figsize=(8,8))
ax  = plt.axes([.1,.1,0.8,0.3], projection=ccrs.PlateCarree())
    
b = [0,0,.0,0]
extent = (lon.min()+b[0],lon.max()+b[1],lat.min()+b[2],lat.max()+b[3])  
    
ax.set_extent( extent )
#ax.coastlines('10m', linewidth=1.,color='gray')
    
    
clevs = np.arange(1,13,.1)
#clevs = np.linspace(temp.quantile(.0),temp.quantile(1)+.5, 100)

arg = {'color':'white', 'lw':1}
tit = 'Urban fraction (-)'

#cmap = plt.get_cmap('coolwarm')
cmap = plt.get_cmap('jet')
    
    
norm = BoundaryNorm(clevs, cmap.N, clip=True)
ax.pcolormesh(lon,lat,temp,norm = norm,cmap= cmap, transform=ccrs.PlateCarree()) #, extend='both')

#ax.contourf(lon,lat,z,levels = clevs,norm = norm,cmap= cmap, transform=ccrs.PlateCarree()) #, extend='both')
#ax.coastlines()

ax.text(1,1.01,'', fontsize=15, fontweight = 'bold',
        ha='right', 
        transform=ax.transAxes)
#add_grid(ax, .25 )
scale_bar(ax, ccrs.PlateCarree(), 100) 
plt.axis('off') 
   
    
    
cbc = 'k'
axc = plt.axes([.1,.18,0.005,0.1])
#axc.set_facecolor("yellow")

cb = mpl.colorbar.ColorbarBase(axc, 
                               cmap=cmap, 
                               norm = norm, 
                               #orientation='horizontal'
                               orientation = 'vertical',
                               format = '%.0f'
                               ) #, extend = 'both',)

cb.ax.tick_params(labelsize=6, color=cbc)

cb.ax.locator_params(nbins=5)
cb.set_label(label='Temp. $\mathrm{^oC}$', fontsize=7)  
#cb.set_label(label=tit, fontsize=10, color=cbc)  
plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=cbc)
plt.setp([cb.ax.get_xticklines(), cb.ax.get_yticklines()], color=cbc)

    
    
