import numpy as np
from datetime import date
from netCDF4 import Dataset
import multiprocessing
from functools import partial
from matplotlib import pyplot as plt
import time

print("Program start")

import storm_functions as storm

filename = 'D:/Robert Docs/College/NIU/GEOG 790 (SP 19)/Project/stormTracking/storm_det_slp'
data = np.load(filename + '.npz', encoding='latin1')
bigListStorms = data['storms']
bigListYear=[]
bigListMonth=[]
bigListDay=[]
bigListHour=[]

dataDir = 'D:/Robert Docs/College/NIU/GEOG 790 (SP 19)/Project/data/'

dataset = 'NARR_MSLET'
			
var = {'NARR_PRMSL': 'prmsl', 
	   'NARR_PRES_SFC': 'pres',
	   'NARR_MSLET': 'mslet'}

# Generate date and hour vectors
yearStart = 1979
yearEnd = 2018 #2018

filename = {'NARR_PRMSL': dataDir + 'prmsl.' + str(yearStart) + '.nc',
            'NARR_PRES_SFC': dataDir + 'pres.sfc.' + str(yearStart) + '.nc',
            'NARR_MSLET': dataDir + 'mslet.' + str(yearStart) + '.nc'}
			
for yr in range(yearStart, yearEnd+1):
	year = np.zeros((0,))
	month = np.zeros((0,))
	day = np.zeros((0,))
	hour = np.zeros((0,))

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

	fileobj.close()
	
	bigListYear.append(year)
	bigListMonth.append(month)
	bigListDay.append(day)
	bigListHour.append(hour)
	
np.savez('storm_det_slp_patched', storms=bigListStorms, year=bigListYear, month=bigListMonth, day=bigListDay, hour=bigListHour)	