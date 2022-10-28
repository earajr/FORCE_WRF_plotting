import numpy as np
import cartopy
from cartopy import crs
from cartopy.feature import NaturalEarthFeature
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import xarray as xr
import os
import ninept_smoother
import calc_gradient
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pint import UnitRegistry

from wrf import (getvar, interplevel, vertcross, CoordPair, ALL_TIMES, to_np, get_cartopy, latlon_coords, cartopy_xlim, cartopy_ylim, extract_times, extract_global_attrs)

# Input WRF out file (this will be given as an argument when running operationally
wrf_fil ="/home/earajr/example_WRF_output/wrfout_d01_2022-05-14_00:00:00"

# Check for existance of WRF out file
if not os.path.exists(wrf_fil):
    raise ValueError("Warning! "+wrf_fil+" does not exist.")

# Read WRF out netcdf
wrf_in= Dataset(wrf_fil)

# Extract the number of times within the WRF file and loop over all times in file
num_times = np.size(extract_times(wrf_in, ALL_TIMES))

for i in np.arange(0, num_times, 1):

# Read 2m temperature, lats and lons

   T2 = getvar(wrf_in, 'T2', timeidx=i)
   lats, lons = latlon_coords(T2)
 
# Calculate T2 gradient

   T2_gradient = calc_gradient.calc_gradient(T2, to_np(lats), to_np(lons))
   T2_gradient_mag_kpm = np.sqrt((T2_gradient[0]**2.0) + (T2_gradient[1]**2.0))
   T2_gradient_mag = T2_gradient_mag_kpm.to("kelvin / kilometer")

## Apply smoothing multiple times to create more user friendly image
#   geop_height_200 = ninept_smoother.smth9(geop_height_200, 0.5, 0.25)
#   geop_height_200 = ninept_smoother.smth9(geop_height_200, 0.5, 0.25)
#   geop_height_200 = ninept_smoother.smth9(geop_height_200, 0.5, 0.25)
#   geop_height_200 = ninept_smoother.smth9(geop_height_200, 0.5, 0.25)

# Read projection from a variable (will be able to detect all possible WRF projections and use them for plotting) 
   cart_proj = get_cartopy(T2)

# Create figure and axes
   fig = plt.figure(figsize=(10,10))
   ax = plt.axes(projection=cart_proj)
   ax.coastlines(linewidth=0.5)

# Plot geopotential height at 10 dam intervals
   T2_gradient_lvl = np.arange(0.0, 0.5, 0.01)
   plt.contourf(lons, lats, T2_gradient_mag, levels=T2_gradient_lvl, cmap='gray_r', transform=crs.PlateCarree())
#  plt.contour(lons, lats, geop_height_200, levels=geop_height_200_lvl, colors='black', transform=crs.PlateCarree())

# Identify whether domain is portrait or landscape

   if np.size(lats[:,0]) < np.size(lats[0,:]):
      portrait = True
   else:
      portrait = False

#  portrait = False

# Create inset colourbar

   if portrait:
      cbbox = inset_axes(ax, '13%', '90%', loc = 7)
      [cbbox.spines[k].set_visible(False) for k in cbbox.spines]
      cbbox.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
      cbbox.set_facecolor([1,1,1,0.7])
#      cbbox.text(0.7,0.5, "200 hPa windspeed (m/s)", rotation=90.0, verticalalignment='center', horizontalalignment='center')
#      cbbox.text(0.85,0.5, "Geopotential height (10 dm spacing)", rotation=90.0, verticalalignment='center', horizontalalignment='center', color='red')
      cbaxes = inset_axes(cbbox, '30%', '95%', loc = 6)
      cb = plt.colorbar(cax=cbaxes, aspect=20)
   else:
      cbbox = inset_axes(ax, '90%', '12%', loc = 8)
      [cbbox.spines[k].set_visible(False) for k in cbbox.spines]
      cbbox.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
      cbbox.set_facecolor([1,1,1,0.7])
#      cbbox.text(0.5,0.3, "200 hPa windspeed (m/s)", verticalalignment='center', horizontalalignment='center')
#      cbbox.text(0.5,0.15, "Geopotential height (10 dm spacing)", verticalalignment='center', horizontalalignment='center', color='red')
      cbaxes = inset_axes(cbbox, '95%', '30%', loc = 9)
      cb = plt.colorbar(cax=cbaxes, orientation='horizontal')

# Add inset timestamp
   tsbox = inset_axes(ax, '95%', '3%', loc = 9)
   [tsbox.spines[k].set_visible(False) for k in tsbox.spines]
   tsbox.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
   tsbox.set_facecolor([1,1,1,1])

   sim_start_time = extract_global_attrs(wrf_in, 'SIMULATION_START_DATE')
   valid_time = str(extract_times(wrf_in, ALL_TIMES)[i])[0:22]

   tsbox.text(0.01, 0.45, "Start date: "+sim_start_time['SIMULATION_START_DATE'], verticalalignment='center', horizontalalignment='left')
   tsbox.text(0.99, 0.45, "Valid_date: "+valid_time, verticalalignment='center', horizontalalignment='right')

## Add wind vectors after thinning.
#   thin = [int(x/15.) for x in lons.shape]
#   ax.quiver(to_np(lons[::thin[0],::thin[1]]), to_np(lats[::thin[0],::thin[1]]), to_np(u_200[::thin[0],::thin[1]]), to_np(v_200[::thin[0],::thin[1]]), pivot='middle', transform=crs.PlateCarree())

# Save image 
   plt.savefig('front_identification.png', bbox_inches='tight')

