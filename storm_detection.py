'''

  Software for the tracking of storms and high-pressure systems

'''

#
# Load required modules
#

import numpy as np
from datetime import date
from netCDF4 import Dataset

from matplotlib import pyplot as plt

import storm_functions as storm

#
# Load in slp data and lat/lon coordinates
#

# Parameters
dataDir = 'D:/Robert Docs/College/NIU/GEOG 790 (SP 19)/Project/data/'

dataset = 'NARR_PRMSL'
			
var = {'NARR_PRMSL': 'prmsl', 
	   'NARR_PRES_SFC': 'pres',
	   'NARR_MSLET': 'mslet'}

# Generate date and hour vectors
yearStart = 1979
yearEnd = 1979 #2018

# Load lat, lon
filename = {'NARR_PRMSL': dataDir + 'prmsl.' + str(yearStart) + '.nc',
            'NARR_PRES_SFC': dataDir + 'pres.sfc.' + str(yearStart) + '.nc',
            'NARR_MSLET': dataDir + 'mslet.' + str(yearStart) + '.nc'}
fileobj = Dataset(filename[dataset], 'r')
lon = fileobj.variables['lon'][:].astype(float)
lat = fileobj.variables['lat'][:].astype(float)
fileobj.close()

print(lon.shape)

# Load slp data
slp = np.zeros((0, lat.shape[0], lat.shape[1]))
year = np.zeros((0,))
month = np.zeros((0,))
day = np.zeros((0,))
hour = np.zeros((0,))
for yr in range(yearStart, yearEnd+1):
	filename = {'NARR_PRMSL': dataDir + 'prmsl.' + str(yr) + '.nc',
				'NARR_PRES_SFC': dataDir + 'pres.sfc.' + str(yr) + '.nc',
				'NARR_MSLET': dataDir + 'mslet.' + str(yr) + '.nc'}
	fileobj = Dataset(filename[dataset], 'r')
	time = fileobj.variables['time'][:]
	time_ordinalDays = time/24. + date(1800,1,1).toordinal()
	year = np.append(year, [date.fromordinal(np.floor(time_ordinalDays[tt]).astype(int)).year for tt in range(len(time))])
	month = np.append(month, [date.fromordinal(np.floor(time_ordinalDays[tt]).astype(int)).month for tt in range(len(time))])
	day = np.append(day, [date.fromordinal(np.floor(time_ordinalDays[tt]).astype(int)).day for tt in range(len(time))])
	hour = np.append(hour, (np.mod(time_ordinalDays, 1)*24).astype(int))
	slp0 = fileobj.variables[var[dataset]][:].astype(float)
	slp = np.append(slp, slp0, axis=0)
	fileobj.close()
	print(yr, slp0.shape[0])

#
# Storm Detection
#

# Initialisation

lon_storms_a = []
lat_storms_a = []
amp_storms_a = []
lon_storms_c = []
lat_storms_c = []
amp_storms_c = []

# Loop over time

T = slp.shape[0]

# Robert: Added this debugging block until I figure out what's going on.
for tS in range(0, 4):
	print(tS)
	#
	# Detect lon and lat coordinates of storms
	#
	lon_storms, lat_storms, amp = storm.detect_storms(slp[tS,:,:], lon, lat, res=2, Npix_min=9, cyc='anticyclonic')
	lon_storms_a.append(lon_storms)
	lat_storms_a.append(lat_storms)
	amp_storms_a.append(amp)
	print(lon_storms)
	print(lat_storms)
	print(amp)
	#
	lon_storms, lat_storms, amp = storm.detect_storms(slp[tS,:,:], lon, lat, res=2, Npix_min=9, cyc='cyclonic')
	lon_storms_c.append(lon_storms)
	lat_storms_c.append(lat_storms)
	amp_storms_c.append(amp)
	#
	# Save as we go
	#
	if tS == 3:
		print('Save data...')
		storms = storm.storms_list(lon_storms_a, lat_storms_a, amp_storms_a, lon_storms_c, lat_storms_c, amp_storms_c)
		np.savez('storm_det_slp', storms=storms, year=year, month=month, day=day, hour=hour)

"""
for tt in range(T):
    #
    print(tt, T)
    #
    # Detect lon and lat coordinates of storms
    #
    lon_storms, lat_storms, amp = storm.detect_storms(slp[tt,:,:], lon, lat, res=2, Npix_min=9, cyc='anticyclonic', globe=False)
    lon_storms_a.append(lon_storms)
    lat_storms_a.append(lat_storms)
    amp_storms_a.append(amp)
    #
    lon_storms, lat_storms, amp = storm.detect_storms(slp[tt,:,:], lon, lat, res=2, Npix_min=9, cyc='cyclonic', globe=False)
    lon_storms_c.append(lon_storms)
    lat_storms_c.append(lat_storms)
    amp_storms_c.append(amp)
    #
    # Save as we go
    #
    if (np.mod(tt, 100) == 0) + (tt == T-1):
        print('Save data...')
    #
    # Combine storm information from all days into a list, and save
    #
        storms = storm.storms_list(lon_storms_a, lat_storms_a, amp_storms_a, lon_storms_c, lat_storms_c, amp_storms_c)
        np.savez('storm_det_slp', storms=storms, year=year, month=month, day=day, hour=hour)
"""
