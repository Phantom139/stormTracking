'''

  Calculate storm census statistics
  for tracked storms

  
  RF NOTE: See http://qingkaikong.blogspot.com/2017/12/use-k-d-tree-to-query-points-part-2-use.html
  
'''

#
# Load required modules
#

import numpy as np
from matplotlib import pyplot as plt
import cartopy
import cartopy.crs as ccrs
from scipy import spatial

### Debug Mode Flag, if set to 1, only does one storm
testMode = False


def storm_in_grid(grid, stormLat, stormLon):
	gridLon_l = grid[0]
	gridLon_u = grid[1]
	gridLat_l = grid[2]
	gridLat_u = grid[3]
	
	return ((stormLat >= gridLat_l and stormLat <= gridLat_u) and (stormLon >= gridLon_l and stormLon <= gridLon_u))
	
#
# Load storm data
#

print("Begin storm census program")

print("Load storm data")
data = np.load('storm_track_slp.npz', encoding='latin1', allow_pickle=True)
storms = data['storms']

# Census grid
plotExtent = [-160, -40, 20, 70]

lonCount = ((plotExtent[1]+360) - (plotExtent[0]+360))
latCount = ((plotExtent[3]) - (plotExtent[2]))

lons = np.linspace(plotExtent[0], plotExtent[1], lonCount)
lats = np.linspace(plotExtent[2], plotExtent[3], latCount)

llon, llat = np.meshgrid(lons, lats)
DIM = llon.shape

if testMode:
	print("Grid:\n- Lons: " + str(plotExtent[0]+360) + " -> " + str(plotExtent[1]+360) + "\n- Lats: " + str(plotExtent[2]) + " -> " + str(plotExtent[3]))
	print("Counts:\n- Lon: " + str(lonCount) + ", Lat: " + str(latCount))
	print("Array Shape: " + str(DIM))
	print("Array:\n" + str(lons) + "\n\n" + str(lats))

tree = spatial.KDTree(list(zip(llon.ravel(), llat.ravel())))

#
# Calculate statistics
#
N = np.zeros(DIM) # Count of track positions
gen = np.zeros(DIM) # Count of genesis positions
term = np.zeros(DIM) # Count of termination positions


if testMode:
	print("Debug mode enabled")
	
	newN = 0
	newGen = 0
	newTerm = 0
	
	storm_lats = storms[0]['lat']
	storm_lons = storms[0]['lon']

	coord_tests = list(zip(storm_lons, storm_lats))		
	distance, index = tree.query(coord_tests)
	closestPoint = tree.data[index]
	
	i, j = np.unravel_index(index, DIM)
	
	for t in range(len(i)):
		storm_lat = storms[0]['lat'][t]
		storm_lon = storms[0]['lon'][t]
		
		print("Point " + str(t) + ": " + str(storm_lon) + ", " + str(storm_lat) + " | Closest: " + str(closestPoint[t]) + " [" + str(i[t]) + ", " + str(j[t]) + "]")
		
		# Check if the storm is inside out bounds
		if(storm_in_grid(plotExtent, storm_lat, storm_lon)):
			N[i[t], j[t]] += 1
			newN += 1
			# Genesis (first location)
			if t == 0:
				gen[i[t], j[t]] += 1
				newGen += 1
			# Termination (last location)
			if t == storms[0]['age']-1:
				term[i[t], j[t]] += 1
				newTerm += 1	

	print("Done processing")
	print("Added: " + str(newN) + " N, " + str(newGen) + " Genesis, " + str(newTerm) + " Termination points")
	print("\n\n")
	print(str(N))
	print(str(gen))
	print(str(term))				
				
else:
	print("Beginning census counting")

	newN = 0
	newGen = 0
	newTerm = 0

	for ed in range(len(storms)):
		storm_lats = storms[ed]['lat']
		storm_lons = storms[ed]['lon']
		
		coord_tests = list(zip(storm_lons, storm_lats))		
		distance, index = tree.query(coord_tests)
		closestPoint = tree.data[index]
		
		i, j = np.unravel_index(index, DIM)
		
		for t in range(len(i)):
			storm_lat = storms[ed]['lat'][t]
			storm_lon = storms[ed]['lon'][t]
			
			# Check if the storm is inside out bounds
			if(storm_in_grid(plotExtent, storm_lat, storm_lon)):
				N[i[t], j[t]] += 1
				newN += 1
				# Genesis (first location)
				if t == 0:
					gen[i[t], j[t]] += 1
					newGen += 1
				# Termination (last location)
				if t == storms[ed]['age']-1:
					term[i[t], j[t]] += 1
					newTerm += 1
						
		if(ed % 200 == 0):
			print("Completed " + str(ed) + " / " + str(len(storms)) + " (" + str((ed/len(storms)) * 100) + "%) -> Added " + str(newN) + " N, " + 
					str(newGen) + " Genesis, " + str(newTerm) + " Termination Points")
			newN = 0
			newGen = 0
			newTerm = 0

	print("Done processing")
	print(str(N))
	print(str(gen))
	print(str(term))

	print("Data Saved")
	np.savez('storm_census', count=N, genesis=gen, termination=term)