'''

  Software for the tracking of storms and high-pressure systems

'''

#
# Load required modules
#

import numpy as np
from datetime import date
from netCDF4 import Dataset
import multiprocessing
from functools import partial
from matplotlib import pyplot as plt
import time

import storm_functions as storm

# Helper function
def run_detection(tS, slp, lon, lat):	
	lon_storms_a = []
	lat_storms_a = []
	amp_storms_a = []
	lon_storms_c = []
	lat_storms_c = []
	amp_storms_c = []	
	# anti-cyclones
	lon_storms, lat_storms, amp = storm.detect_storms(slp[tS,:,:], lon, lat, res=2, Npix_min=9, cyc='anticyclonic')
	lon_storms_a.append(lon_storms)
	lat_storms_a.append(lat_storms)
	amp_storms_a.append(amp)
	# cyclones
	lon_storms, lat_storms, amp = storm.detect_storms(slp[tS,:,:], lon, lat, res=2, Npix_min=9, cyc='cyclonic')
	lon_storms_c.append(lon_storms)
	lat_storms_c.append(lat_storms)
	amp_storms_c.append(amp)
	# Write out
	storms = storm.storms_list(lon_storms_a, lat_storms_a, amp_storms_a, lon_storms_c, lat_storms_c, amp_storms_c)
	return storms

#
# Load in slp data and lat/lon coordinates
#

# Parameters
## NOTE: MAKE SURE YOU EDIT THIS LINE!!!!
## THIS IS WHERE THE PROGRAM WILL LOOK FOR DATA
dataDir = 'D:/Robert Docs/College/NIU/GEOG 790 (SP 19)/Project/data/'

dataset = 'NARR_MSLET'
			
var = {'NARR_PRMSL': 'prmsl', 
	   'NARR_PRES_SFC': 'pres',
	   'NARR_MSLET': 'mslet'}

# Generate date and hour vectors
yearStart = 1979
yearEnd = 2018 #2018

# Load lat, lon
filename = {'NARR_PRMSL': dataDir + 'prmsl.' + str(yearStart) + '.nc',
            'NARR_PRES_SFC': dataDir + 'pres.sfc.' + str(yearStart) + '.nc',
            'NARR_MSLET': dataDir + 'mslet.' + str(yearStart) + '.nc'}
fileobj = Dataset(filename[dataset], 'r')
lon = fileobj.variables['lon'][:].astype(float)
lat = fileobj.variables['lat'][:].astype(float)
fileobj.close()

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
	timeAr = fileobj.variables['time'][:]
	time_ordinalDays = timeAr/24. + date(1800,1,1).toordinal()
	year = np.append(year, [date.fromordinal(np.floor(time_ordinalDays[tt]).astype(int)).year for tt in range(len(timeAr))])
	month = np.append(month, [date.fromordinal(np.floor(time_ordinalDays[tt]).astype(int)).month for tt in range(len(timeAr))])
	day = np.append(day, [date.fromordinal(np.floor(time_ordinalDays[tt]).astype(int)).day for tt in range(len(timeAr))])
	hour = np.append(hour, (np.mod(time_ordinalDays, 1)*24).astype(int))
	slp0 = fileobj.variables[var[dataset]][:].astype(float)
	slp = np.append(slp, slp0, axis=0)
	fileobj.close()
	print(yr)

#
# Storm Detection
#


# Loop over time
bigListStorms = []
T = slp.shape[0]
totalTime = 0
print("Size of T: " + str(T))
processStart = time.time()
for tS in range(500): #T):
	if(tS == 0):
		print("Processing " + str(tS+1) + "/" + str(T+1))
	else:
		avg = totalTime / (tS+1)
		est = avg * (T - tS)
		print("Processing " + str(tS+1) + "/" + str(T+1) + ": Last Step: " + '{0:.2f}'.format(total) + "s, Avg. Time: " + '{0:.2f}'.format(avg) + "s, Est. Time Left: " + time.strftime("%H:%M:%S", time.gmtime(est)))
		
	timeStart = time.time()
	storms = run_detection(tS, slp, lon, lat)
	bigListStorms.append(storms)
	timeEnd = time.time()
	
	total = timeEnd - timeStart
	totalTime += total
	
	if(tS % 50 == 0 or tS == T-1):
		np.savez('storm_det_slp', storms=bigListStorms, year=year, month=month, day=day, hour=hour)
processEnd = time.time()
elapsed = processEnd - processStart

print("Program completed. Total Elapsed Time: " + time.strftime("%H:%M:%S", time.gmtime(elapsed)))