'''

  Software for the tracking of storms
  based on detected storm position data.

'''

# Load required modules

import numpy as np
import storm_functions as storm

#
# Automated storm tracking
#

# Load in detected positions and date/hour information
filename = 'D:/Robert Docs/College/NIU/GEOG 790 (SP 19)/Project/stormTracking/storm_det_slp'
data = np.load(filename + '.npz', encoding='latin1')
det_storms = data['storms']
year = data['year']
month = data['month']
day = data['day']
hour = data['hour']

print(det_storms.shape)
print(year.shape)

# Initialize storms discovered at first time step

storms = storm.storms_init(det_storms, year, month, day, hour)

# Stitch storm tracks together at future time steps

T = len(det_storms) # number of time steps
for tt in range(0, T):
	# Track storms from time step tt-1 to tt and update corresponding tracks and/or create new storms
	storms = storm.track_storms(storms, det_storms, tt, year, month, day, hour, dt=3)

# Add keys for storm age and flag if storm was still in existence at end of run
for ed in range(len(storms)):
    storms[ed]['age'] = len(storms[ed]['lon'])

# Robert: Added Classification here
storms = storm.classify_storms(storms)	
	
# Strip storms based on track lengths (dt = 3hr, d_tot_min = 500km, dur_min = 24hr)
storms = storm.strip_storms(storms, dt=3, d_tot_min=500., d_ratio=0.6, dur_min=24)

# Save tracked storm data
np.savez('storm_track_slp', storms=storms)