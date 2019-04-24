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
			
for yr in data['year']:
	bigListYear.extend(yr)
for yr in data['month']:
	bigListMonth.extend(yr)
for yr in data['day']:
	bigListDay.extend(yr)
for yr in data['hour']:
	bigListHour.extend(yr)	
	
print(bigListYear)	

np.savez('storm_det_slp_patched', storms=bigListStorms, year=bigListYear, month=bigListMonth, day=bigListDay, hour=bigListHour)	