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
import os
import storm_functions as storm
import itertools
import time

# Helper function
def run_detection(slp, idx, size, lon, lat):	
	processStart = time.time()
	lon_storms_a = []
	lat_storms_a = []
	amp_storms_a = []
	lon_storms_c = []
	lat_storms_c = []
	amp_storms_c = []	
	# anti-cyclones
	lon_storms, lat_storms, amp = storm.detect_storms(slp, lon, lat, res=2, Npix_min=9, cyc='anticyclonic')
	lon_storms_a.append(lon_storms)
	lat_storms_a.append(lat_storms)
	amp_storms_a.append(amp)
	# cyclones
	lon_storms, lat_storms, amp = storm.detect_storms(slp, lon, lat, res=2, Npix_min=9, cyc='cyclonic')
	lon_storms_c.append(lon_storms)
	lat_storms_c.append(lat_storms)
	amp_storms_c.append(amp)
	# Write out
	storms = storm.storms_list(lon_storms_a, lat_storms_a, amp_storms_a, lon_storms_c, lat_storms_c, amp_storms_c)
	processEnd = time.time()
	print("Time step completed (" + str(idx) + " / " + str(size) + "), Elapsed Time: " + time.strftime("%H:%M:%S", time.gmtime(processEnd - processStart)))
	return {idx: storms}

#
# Load in slp data and lat/lon coordinates
#
print("Program Start...")

# Parameters
## NOTE: MAKE SURE YOU EDIT THIS LINE!!!!
## THIS IS WHERE THE PROGRAM WILL LOOK FOR DATA
dataDir = '/media/robert/HDD2/790Project/Data/' #'D:/Robert Docs/College/NIU/GEOG 790 (SP 19)/Project/data/'

dataset = 'NARR_MSLET'
			
var = {'NARR_PRMSL': 'prmsl', 
	   'NARR_PRES_SFC': 'pres',
	   'NARR_MSLET': 'mslet'}

# Generate date and hour vectors
yearStart = 1979 # 1979
yearEnd = 2020 #2018

# Load lat, lon
filename = {'NARR_PRMSL': dataDir + 'prmsl.' + str(yearStart) + '.nc',
            'NARR_PRES_SFC': dataDir + 'pres.sfc.' + str(yearStart) + '.nc',
            'NARR_MSLET': dataDir + 'mslet.' + str(yearStart) + '.nc'}
			
print("Loading in first netCDF file to populate arrays... " + str(filename[dataset]))
fileobj = Dataset(filename[dataset], 'r')
lon = fileobj.variables['lon'][:].astype(float)
lat = fileobj.variables['lat'][:].astype(float)
fileobj.close()

bigListStorms = []

print("Spawning a multiprocessing pool, " + str(os.cpu_count()) + " processors detected.")
pool = multiprocessing.Pool(os.cpu_count())

print("Entering loop, beginning timer.")
fullStart = time.time()
yListTime = []

# Create empty arrays to hold the data
year = np.zeros((0,))
month = np.zeros((0,))
day = np.zeros((0,))
hour = np.zeros((0,))

for yr in range(yearStart, yearEnd+1):
	save_storms = {}
	current_list_storms = []

	innerTimeStart = time.time()

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
	fileobj.close()
	
	# Run the multiprocessed detection function, merge the returned dictionary object with our local one.
	out_dict = pool.starmap(run_detection, zip(slp0, np.arange(slp0.shape[0]), itertools.repeat(slp0.shape[0]), itertools.repeat(lon), itertools.repeat(lat)))
	ro_inner = time.time()
	for iDict in out_dict: #Note: starmap returns as a 0-length list ie [return] with the Nth element being each input
		save_storms = {**save_storms, **iDict} 
	# Order the final list correctly
	for i in range(slp0.shape[0]):
		current_list_storms.append(save_storms[i])
	ro_outer = time.time() - ro_inner
	print(str(yr) + " List Reordering Completed, Elapsed Time: " + time.strftime("%H:%M:%S", time.gmtime(ro_outer)))
	# Append to the big list, then save.
	bigListStorms.append(current_list_storms)
	innerTimeEnd = time.time()
	yListTime.append((innerTimeEnd - innerTimeStart))
	print(str(yr) + " Completed, Elapsed Time: " + time.strftime("%H:%M:%S", time.gmtime(innerTimeEnd - innerTimeStart)))
	
	np.savez('storm_det_slp', storms=bigListStorms, year=year, month=month, day=day, hour=hour)

pool.close() 
pool.join()

fullEnd = time.time()
print("Program Completed, Elapsed Time: " + time.strftime("%H:%M:%S", time.gmtime(fullEnd - fullStart)))

for i, yi in enumerate(yListTime):
	print("Processing Time (" + str(yearStart + i) + "): " + time.strftime("%H:%M:%S", time.gmtime(yi)))
