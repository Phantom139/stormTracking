import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import cartopy
import cartopy.crs as ccrs
import scipy
from scipy import spatial
from scipy import stats
import storm_functions as storm

# Load storm data
det_storms = np.load('storm_det_slp.npz', allow_pickle=True)
stormPosn = det_storms['storms']

tracked_storms = np.load('storm_track_slp.npz', allow_pickle=True)
storms = tracked_storms['storms']

print("Debug Script\n\n")
print("det_storms\n\n")
print(stormPosn)
print("\n\ntracked_storms\n\n")
print(storms)

def calculate_bergeron_test(stormObject):
	bergeronList = []
	for i in range(len(stormObject['amp'])):
		pressures = stormObject['amp'][i:i+8]
		latitudes = stormObject['lat'][i:i+8]
		
		angFactor = np.sin(np.radians(60)) / np.sin(np.radians(latitudes))
		
		diffs = np.diff(pressures) / 100 #Convert Pa to mb
		negatives = np.sum(val for val in diffs if val < 0)
		
		
		bFactor = (-1 * np.sum(negatives)) / (24)
		
		final = bFactor * angFactor
		
		bergeronList.append(bFactor)
		
	if(any(b >= 1 for b in bergeronList)):	
		print("Bergeron > 1 (" + str(bergeronList) + ") " + str(stormObject['amp']))
	
	
for ed in range(len(storms)):
	calculate_bergeron_test(storms[ed])