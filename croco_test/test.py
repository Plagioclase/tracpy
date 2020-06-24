#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 15:37:38 2020

@author: noam
"""

import numpy as np
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
from tracpy.tracpy_class import Tracpy
import matplotlib.pyplot as plt
import matplotlib
import datetime
from mpl_toolkits.basemap import Basemap

# File location
loc = '../../croco_sim/Run_WIOMonClim_Coburg10/CROCO_FILES/croco_12h.nc'

# Number of days to run the drifters.
ndays = 14

# Start date in date time formatting
date = datetime.datetime(2005, 12, 1, 0)

# Time between outputs
tseas = 12*3600  # 4 hours between outputs, in seconds

# Time units
time_units = 'seconds since 2000-01-01 00:00:00'

# Sets a smaller limit than between model outputs for when to force
# interpolation if hasn't already occurred.
nsteps = 5
# Controls the sampling frequency of the drifter tracks.
N = 4

# Use ff = 1 for forward in time and ff = -1 for backward in time.
ff = 1

ah = 10.  # m^2/s
av = 0.  # m^2/s

# turbulence/diffusion flag
doturb = 1

# simulation name, used for saving results into netcdf file
name = 'temp'

# for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
do3d = 1

# Choose method for vertical placement of drifters
z0 = np.linspace(-5,-5,10000)#-3000*np.ones(100)#-5*np.ones(900) # I know the size from checking #'s after eliminating those outside
#           domain ' #'z' #'salt' #'s' 
num_layers = 50
zpar = 'fromZeta' # 29 #-10 #grid['km']-1 # 30 #grid['km']-1

# #### 3D Sample Options ####
# # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
# do3d = 1

# ## Choose method for vertical placement of drifters
# z0 = np.zeros(676) # I know the size from checking #'s after eliminating those outside domain ' #'z' #'salt' #'s' 
# num_layers = 30
# zpar = 'fromZeta' #num_layers-1 # 29 #-10 #grid['km']-1 # 30 #grid['km']-1
# ####

proj = tracpy.tools.make_proj('WIO')

# file, tinds = tracpy.inout.setupROMSfiles(loc, date, ff, 30, time_units)
roms_grid = tracpy.inout.readgrid(loc, proj)
tp = Tracpy(loc, roms_grid, name=name, tseas=tseas, ndays=ndays, nsteps=nsteps,
            N=N, ff=ff, ah=ah, av=av, doturb=doturb, do3d=do3d, z0=z0, zpar=zpar, time_units=time_units)

# lon0, lat0 = np.meshgrid(np.linspace(49.25,50.75,50), \
#                             np.linspace(-15.0,-12.5,50)) # whole domain, 20 km
lon0 = np.array([55.36, 55.80, 46.18, 50.59, 57.69, 41.00, 72.4, 59.61, 45.14, 53.33])
lat0 = np.array([-4.56, -4.24, -9.33, -15.38, -19.96, -14.65, -7.48, -16.63, -13.04, -5.43])

# Eliminate points that are outside domain or in masked areas
lon0, lat0 = tracpy.tools.check_points(lon0, lat0, tp.grid)

lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)

fig = plt.figure(figsize=(20, 20), dpi=100)
fig, ax = tracpy.plotting.background(roms_grid, fig=fig, extent=[35, 60, -20, 0],
               col='lightgrey', halpha=1, outline=[1, 1, 0, 1], res='50m')
tracpy.plotting.tracks(lonp, latp, tp.name, roms_grid, fig=fig, ax=ax)
