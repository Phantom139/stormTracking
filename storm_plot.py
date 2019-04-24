'''

  Plot storm tracks

'''

# Load required modules

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import cartopy
import cartopy.crs as ccrs

# Load storm data
det_storms = np.load('storm_det_slp.npz')
stormPosn = det_storms['storms']

tracked_storms = np.load('storm_track_slp.npz')
storms = tracked_storms['storms']

# Settings
plotExtent = [-125, -70, 23, 49] # Where to show the plot (LonMin, LonMax, LatMin, LatMax)

def get_projection_object():
    # Cartopy has a "globe" object to define more projection standards, we create this first
    globe = ccrs.Globe(ellipse=None,
                       semimajor_axis=6370000,
                       semiminor_axis=6370000,
                       nadgrids="@null")    
    # Now we can create the projection object
    cen_lat  = float(50.0)
    cen_lon  = float(-107.0)
    std_pll  = [50.0]
    cutoff   = -30.0
    projObj = ccrs.LambertConformal(central_latitude=cen_lat, 
                                    central_longitude=cen_lon,
                                    standard_parallels=std_pll,
                                    globe=globe,
                                    cutoff=cutoff)
    
    return projObj

def plot_positions(stormObj):
	# Grab the projection object from the netCDF file
	projObj = get_projection_object()    
	# Create our figure and axis object
	fig = plt.figure(figsize=(12,9))
	ax = plt.axes(projection=projObj)

	ax.set_extent([-125, -70, 23, 60], crs=ccrs.PlateCarree())   

	# Draw our plot, coastlines first, then the contours.
	states = cartopy.feature.NaturalEarthFeature(category='cultural',
												name='admin_1_states_provinces_lakes',
												scale='110m',
												facecolor='none')
	ax.add_feature(states, edgecolor='k')    
	ax.coastlines()

	unassigned = list(range(stormObj[0][0]['N']))
	for i in unassigned:
		lat, lon = stormObj[0][0]['lat'], stormObj[0][0]['lon']
		plt.plot(lon, lat, 'o', linewidth=1, alpha=0.35, transform=ccrs.PlateCarree())
	unassigned2 = list(range(stormObj[1][0]['N']))
	for i in unassigned2:
		lat, lon = stormObj[1][0]['lat'], stormObj[1][0]['lon']
		plt.plot(lon, lat, 'ro', linewidth=1, alpha=0.35, transform=ccrs.PlateCarree())		
	# Show the plot.
	plt.title('Storm Points')
	plt.savefig('figures/storm_points', bbox_inches='tight', pad_inches=0.05, dpi=300)			
	
def plot_tracks(stormObj):
	# Grab the projection object from the netCDF file
	projObj = get_projection_object()    
	# Create our figure and axis object
	fig = plt.figure(figsize=(12,9))
	ax = plt.axes(projection=projObj)
	#ax = plt.axes(projection=ccrs.PlateCarree())
	ax.set_extent([-125, -70, 23, 60], crs=ccrs.PlateCarree())   

	# Draw our plot, coastlines first, then the contours.
	states = cartopy.feature.NaturalEarthFeature(category='cultural',
												name='admin_1_states_provinces_lakes',
												scale='110m',
												facecolor='none')
	ax.add_feature(states, edgecolor='k')    
	ax.coastlines()

	for ed in range(len(stormObj)):
		if (stormObj[ed]['type'] == 'cyclonic'):
			lon, lat = stormObj[ed]['lon'], stormObj[ed]['lat']
			
			lon[lon <  0] = 360 + lon[lon < 0]
			
			plt.plot(lon, lat, 'r-', linewidth=1, alpha=0.35, transform=ccrs.PlateCarree())
			plt.plot(lon[0], lat[0], 'bo', alpha=0.5, markeredgewidth=0, transform=ccrs.PlateCarree())
			plt.plot(lon[-1], lat[-1], 'ro', alpha=0.5, markeredgewidth=0, transform=ccrs.PlateCarree())
		
	# Show the plot.
	plt.plot(0, 0, 'bo', alpha=1, markeredgewidth=0, transform=ccrs.PlateCarree(), label="Cyclogenesis")
	plt.plot(0, 0, 'ro', alpha=1, markeredgewidth=0, transform=ccrs.PlateCarree(), label="Cyclosis")
	plt.legend()
	plt.title('Storm Tracks')
	plt.savefig('figures/storm_tracks', bbox_inches='tight', pad_inches=0.05, dpi=300)
    
def plot_frequency(stormObj, months='A', years='A', cyclones_to_include = 'A'):
	monthList = []
	yearList = []
	cycloneList = []
	
	if months == 'A':
		monthList = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
	else:
		monthList = months
	
	if years == 'A':
		yearList = np.arange(1979, 2019, 1)
	else:
		yearList = years
		
	if cyclones_to_include == 'A':
		cycloneList = ['Clipper', 'Colorado', 'Other']
	else:
		cycloneList = cyclones_to_include
		
	lonRange = np.linspace(plotExtent[0], plotExtent[1], 200)
	latRange = np.linspace(plotExtent[2], plotExtent[3], 200)
	lons, lats = np.meshgrid(lonRange, latRange)
		
	"""
	Robert Notes for Group:
		- The above creates a template listable call where you can pass a list of months, years, and cyclone types
		   to be saved in a list instance. As we're dealing with loop iteration we can take advantage of python's
		   list membership capability to make testing easy, here's an example:
		   
	# Make these first two lines a "standard" practice.
	for ed in range(len(stormObj)):
		if (stormObj[ed]['type'] == 'cyclonic'):
			# Test months
			if(storms[ed]['month'][0] in monthList):
				# Test years
				if(storms[ed]['year'][0] in yearList):
					# Test classification
					if(storms[ed]['classification'] in cycloneList):
						# Anything inside this if statement will have met all three of the above conditions.
	"""
	
# List of function calls, uncomment out the plots you want
#plot_positions(stormPosn)
plot_tracks(storms)

"""
# Plot storm tracks

# Example looking down over the North Pole
plt.figure()
plt.clf()
#
plt.subplot(1,2,1)
proj = bm.Basemap(projection='npstere',boundinglat=25,lon_0=300,resolution='l')
proj.drawcoastlines(linewidth=0.5)
for ed in range(len(storms)):
    # Select for: cyclonic storms which exist solely during November-March
    if (storms[ed]['type'] == 'cyclonic') * ((storms[ed]['month'][0] >= 11)+(storms[ed]['month'][0] <= 3)) * ((storms[ed]['month'][-1]>=11)+(storms[ed]['month'][-1] <= 3)) * (storms[ed]['year'][0] >= 2000) * (storms[ed]['year'][-1] <= 2010):
        lonproj, latproj = proj(storms[ed]['lon'], storms[ed]['lat'])
        plt.plot(lonproj, latproj, 'r-', linewidth=1, alpha=0.25)
plt.title('Storm tracks (Nov-Mar, 2000-2010)')
#
plt.subplot(1,2,2)
proj = bm.Basemap(projection='npstere',boundinglat=25,lon_0=300,resolution='l')
proj.drawcoastlines(linewidth=0.5)
for ed in range(len(storms)):
    # Select for: cyclonic storms which exist solely during November-March
    if (storms[ed]['type'] == 'cyclonic') * ((storms[ed]['month'][0] >= 11)+(storms[ed]['month'][0] <= 3)) * ((storms[ed]['month'][-1]>=11)+(storms[ed]['month'][-1] <= 3)) * (storms[ed]['year'][0] >= 1950) * (storms[ed]['year'][-1] <= 1960):
        lonproj, latproj = proj(storms[ed]['lon'], storms[ed]['lat'])
        plt.plot(lonproj, latproj, 'r-', linewidth=1, alpha=0.25)
plt.title('Storm tracks (Nov-Mar, 1950-1960)')
plt.savefig('figures/storm_tracks_NorthernHemisphere', bbox_inches='tight', pad_inches=0.05, dpi=300)

# South Pole
plt.clf()
#
plt.subplot(1,2,1)
proj = bm.Basemap(projection='spstere',boundinglat=-25,lon_0=300,resolution='l')
proj.drawcoastlines(linewidth=0.5)
for ed in range(len(storms)):
    # Select for: cyclonic storms which exist solely during November-March
    if (storms[ed]['type'] == 'cyclonic') * ((storms[ed]['month'][0] >= 11)+(storms[ed]['month'][0] <= 3)) * ((storms[ed]['month'][-1]>=11)+(storms[ed]['month'][-1] <= 3)) * (storms[ed]['year'][0] >= 2000) * (storms[ed]['year'][-1] <= 2010):
        lonproj, latproj = proj(storms[ed]['lon'], storms[ed]['lat'])
        plt.plot(lonproj, latproj, 'r-', linewidth=1, alpha=0.25)
plt.title('Storm tracks (Nov-Mar, 2000-2010)')
#
plt.subplot(1,2,2)
proj = bm.Basemap(projection='spstere',boundinglat=-25,lon_0=300,resolution='l')
proj.drawcoastlines(linewidth=0.5)
for ed in range(len(storms)):
    # Select for: cyclonic storms which exist solely during November-March
    if (storms[ed]['type'] == 'cyclonic') * ((storms[ed]['month'][0] >= 11)+(storms[ed]['month'][0] <= 3)) * ((storms[ed]['month'][-1]>=11)+(storms[ed]['month'][-1] <= 3)) * (storms[ed]['year'][0] >= 1950) * (storms[ed]['year'][-1] <= 1960):
        lonproj, latproj = proj(storms[ed]['lon'], storms[ed]['lat'])
        plt.plot(lonproj, latproj, 'r-', linewidth=1, alpha=0.25)
plt.title('Storm tracks (Nov-Mar, 1950-1960)')
plt.savefig('figures/storm_tracks_SouthernHemisphere', bbox_inches='tight', pad_inches=0.05, dpi=300)

# Regional zoom in on the NW Atlantic
plt.figure()
plt.clf()
proj = bm.Basemap(llcrnrlon=-100.,llcrnrlat=15.,urcrnrlon=0.,urcrnrlat=65., projection='lcc',lat_1=20.,lat_2=40.,lon_0=-60., resolution ='l',area_thresh=1000.)
proj.drawcoastlines(linewidth=0.5)
for ed in range(len(storms)):
    # Select for: cyclonic storms which exist solely during November-March
    if (storms[ed]['type'] == 'cyclonic') * ((storms[ed]['month'][0] >= 11)+(storms[ed]['month'][0] <= 3)) * ((storms[ed]['month'][-1]>=11)+(storms[ed]['month'][-1] <= 3)) * (storms[ed]['year'][0] >= 2000) * (storms[ed]['year'][-1] <= 2005):
        lonproj, latproj = proj(storms[ed]['lon'], storms[ed]['lat'])
        plt.plot(lonproj, latproj, '-', linewidth=1, alpha=0.6)
plt.title('Storm tracks (2000-2005)')
plt.savefig('figures/storm_tracks_NorthwestAtlantic', bbox_inches='tight', pad_inches=0.05, dpi=300)

# Genesis and termination locations
plt.clf()
proj = bm.Basemap(llcrnrlon=-100.,llcrnrlat=15.,urcrnrlon=0.,urcrnrlat=65., projection='lcc',lat_1=20.,lat_2=40.,lon_0=-60., resolution ='l',area_thresh=1000.)
proj.drawcoastlines(linewidth=0.5)
for ed in range(len(storms)):
    # Select for: cyclonic storms which exist solely during November-March
    if (storms[ed]['type'] == 'cyclonic') * ((storms[ed]['month'][0] >= 11)+(storms[ed]['month'][0] <= 3)) * ((storms[ed]['month'][-1]>=11)+(storms[ed]['month'][-1] <= 3)) * (storms[ed]['year'][0] >= 2000) * (storms[ed]['year'][-1] <= 2010):
        lonproj, latproj = proj(storms[ed]['lon'], storms[ed]['lat'])
        plt.plot(lonproj[0], latproj[0], 'bo', alpha=0.5, markeredgewidth=0) # Genesis point
        plt.plot(lonproj[-1], latproj[-1], 'ro', alpha=0.5, markeredgewidth=0) # End point
plt.legend(['Genesis', 'Termination'], loc='lower right')
plt.title('Storm genesis and termination locations (2000-2010)')
plt.savefig('figures/storm_genesisAndTermination_NorthwestAtlantic', bbox_inches='tight', pad_inches=0.05, dpi=300)

# Only storms which originated in the Mar 1855 or Mar 2010, with central pressures indicated
plt.figure()
plt.clf()
plt.subplot(1,2,1)
proj = bm.Basemap(llcrnrlon=-80.,llcrnrlat=30.,urcrnrlon=-20.,urcrnrlat=65., projection='lcc',lat_1=20.,lat_2=40.,lon_0=-60., resolution ='l',area_thresh=1000.)
proj.drawcoastlines(linewidth=0.5)
for ed in range(len(storms)):
    # Select for: cyclonic storms which started solely during Mar of 1855
    if (storms[ed]['type'] == 'cyclonic') * (storms[ed]['month'][0] == 3)*(storms[ed]['year'][0] == 1855):
        lonproj, latproj = proj(storms[ed]['lon'], storms[ed]['lat'])
        plt.plot(lonproj, latproj, '-', linewidth=1, zorder=10)
        plt.scatter(lonproj, latproj, c=storms[ed]['amp'], zorder=20)
H = plt.colorbar()
H.set_label('Central pressure [Pa]')
plt.title('Storms starting in Mar''1855')
plt.subplot(1,2,2)
proj = bm.Basemap(llcrnrlon=-80.,llcrnrlat=30.,urcrnrlon=-20.,urcrnrlat=65., projection='lcc',lat_1=20.,lat_2=40.,lon_0=-60., resolution ='l',area_thresh=1000.)
proj.drawcoastlines(linewidth=0.5)
for ed in range(len(storms)):
    # Select for: cyclonic storms which started solely during Mar of 2010
    if (storms[ed]['type'] == 'cyclonic') * (storms[ed]['month'][0] == 3)*(storms[ed]['year'][0] == 2010):
        lonproj, latproj = proj(storms[ed]['lon'], storms[ed]['lat'])
        plt.plot(lonproj, latproj, '-', linewidth=1, zorder=10)
        plt.scatter(lonproj, latproj, c=storms[ed]['amp'], zorder=20)
H = plt.colorbar()
H.set_label('Central pressure [Pa]')
plt.title('Storms starting in Mar''2010')
# plt.savefig('figures/storm_CentralPressures_March_1855_2010', bbox_inches='tight', pad_inches=0.05, dpi=300)
"""