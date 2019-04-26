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
import storm_functions as storm

# Load storm data
det_storms = np.load('storm_det_slp.npz')
stormPosn = det_storms['storms']

tracked_storms = np.load('storm_track_slp.npz')
storms = tracked_storms['storms']

# Settings
plotExtent = [-125, -70, 23, 60] # Where to show the plot (LonMin, LonMax, LatMin, LatMax)

zoneColors = {
"Clipper": "#3687FF",
"Northwest": "#6C179E",
"Colorado": "#F3A200",
"GreatBasin": "#A45E45",
"GulfOfMexico": "#4B9F0B",
"EastCoast": "#DE55D9",
"Other": "#666666",
}
typeStr = {
"Clipper": "Alberta Clipper",
"Northwest": "Northwestern",
"Colorado": "Colorado Low",
"GreatBasin": "Great Basin Low",
"GulfOfMexico": "Gulf of Mexico Low",
"EastCoast": "East Coast Low",
"Other": "Other Lows",
}

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
	
def plot_zones():
	clipper_zone = [[-115,50],
					[-105,50],
					[-105,55],
					[-110,55],
					[-110,60],
					[-125,60],
					[-125,55],
					[-115,55]]
	northwst_zone = [[-125, 60],
					[-115,60],
					[-115,65],
					[-125,65]]
	colorado_zone = [[-105, 40],
					[-100,40],
					[-100,35],
					[-105,35]]
	basin_zone = 	[[-120, 45],
					[-115,45],
					[-115,35],
					[-120,35]]
	gom_zone =      [[-100, 25],
					[-90,25],
					[-90,30],
					[-100,30]]
	eastcst_zone =  [[-80, 30],
					[-75,30],
					[-75,35],
					[-65,35],
					[-65,45],
					[-70,45],
					[-70,40],
					[-80,40]]
					
	projObj = get_projection_object()    
	# Create our figure and axis object
	fig = plt.figure(figsize=(12,9))
	ax = plt.axes(projection=projObj)
	#ax = plt.axes(projection=ccrs.PlateCarree())
	ax.set_extent(plotExtent, crs=ccrs.PlateCarree())  
	#ax.set_xlim([-125, -70])
	#ax.set_ylim([23, 49]) 	

	# Draw our plot, coastlines first, then the contours.
	states = cartopy.feature.NaturalEarthFeature(category='cultural',
												name='admin_1_states_provinces_lakes',
												scale='110m',
												facecolor='none')
	ax.add_feature(states, edgecolor='k')    
	ax.coastlines()					
					
	clipperPoly = mplPath.Path(clipper_zone)
	northwestPoly = mplPath.Path(northwst_zone)
	coloradoPoly = mplPath.Path(colorado_zone)
	basinPoly = mplPath.Path(basin_zone)
	gomPoly = mplPath.Path(gom_zone)
	eastCoastPoly = mplPath.Path(eastcst_zone)	
	
	patch1 = matplotlib.patches.PathPatch(clipperPoly, facecolor=zoneColors["Clipper"], transform=ccrs.PlateCarree(), label="Alberta Clipper")
	patch2 = matplotlib.patches.PathPatch(northwestPoly, facecolor=zoneColors["Northwest"], transform=ccrs.PlateCarree(), label="Northwestern")
	patch3 = matplotlib.patches.PathPatch(coloradoPoly, facecolor=zoneColors["Colorado"], transform=ccrs.PlateCarree(), label="Colorado Low")
	patch4 = matplotlib.patches.PathPatch(basinPoly, facecolor=zoneColors["GreatBasin"], transform=ccrs.PlateCarree(), label="Great Basin Low")
	patch5 = matplotlib.patches.PathPatch(gomPoly, facecolor=zoneColors["GulfOfMexico"], transform=ccrs.PlateCarree(), label="Gulf of Mexico Low")
	patch6 = matplotlib.patches.PathPatch(eastCoastPoly, facecolor=zoneColors["EastCoast"], transform=ccrs.PlateCarree(), label="East Coast Low")
	ax.add_patch(patch1)	
	ax.add_patch(patch2)
	ax.add_patch(patch3)
	ax.add_patch(patch4)
	ax.add_patch(patch5)
	ax.add_patch(patch6)
	
	plt.legend()
	plt.title('Storm Classification Zones')
	plt.savefig('figures/storm_classes', bbox_inches='tight', pad_inches=0.05, dpi=300)	
	
def plot_tracks(stormObj):
	# Grab the projection object from the netCDF file
	projObj = get_projection_object()    
	# Create our figure and axis object
	fig = plt.figure(figsize=(12,9))
	ax = plt.axes(projection=projObj)
	#ax = plt.axes(projection=ccrs.PlateCarree())
	ax.set_extent([-125, -70, 23, 70], crs=ccrs.PlateCarree())   

	# Draw our plot, coastlines first, then the contours.
	states = cartopy.feature.NaturalEarthFeature(category='cultural',
												name='admin_1_states_provinces_lakes',
												scale='110m',
												facecolor='none')
	ax.add_feature(states, edgecolor='k')    
	ax.coastlines()
	for ed in range(len(stormObj)):
		if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == 'Clipper'):
			if(stormObj[ed]['month'][0] in [1, 2, 3, 4, 10, 11, 12]):
				lon, lat = stormObj[ed]['lon'], stormObj[ed]['lat']
				
				lon[lon <  0] = 360 + lon[lon < 0]
				
				plt.plot(lon, lat, 'r-', linewidth=0.5, alpha=0.35, transform=ccrs.PlateCarree())
				plt.plot(lon[0], lat[0], 'bo', alpha=0.25, markeredgewidth=0, transform=ccrs.PlateCarree())
				plt.plot(lon[-1], lat[-1], 'ro', alpha=0.25, markeredgewidth=0, transform=ccrs.PlateCarree())
	# Show the plot.
	plt.plot(0, 0, 'bo', alpha=1, markeredgewidth=0, transform=ccrs.PlateCarree(), label="Cyclogenesis")
	plt.plot(0, 0, 'ro', alpha=1, markeredgewidth=0, transform=ccrs.PlateCarree(), label="Cyclosis")
	plt.legend()
	plt.title('Storm Tracks')
	plt.savefig('figures/storm_tracks', bbox_inches='tight', pad_inches=0.05, dpi=300)
	
def plot_bomb_tracks(stormObj):
	# Grab the projection object from the netCDF file
	projObj = get_projection_object()    
	# Create our figure and axis object
	fig = plt.figure(figsize=(12,9))
	ax = plt.axes(projection=projObj)
	#ax = plt.axes(projection=ccrs.PlateCarree())
	ax.set_extent([-125, -70, 23, 70], crs=ccrs.PlateCarree())   

	# Draw our plot, coastlines first, then the contours.
	states = cartopy.feature.NaturalEarthFeature(category='cultural',
												name='admin_1_states_provinces_lakes',
												scale='110m',
												facecolor='none')
	ax.add_feature(states, edgecolor='k')    
	ax.coastlines()
	for ed in range(len(stormObj)):
		if (stormObj[ed]['type'] == 'cyclonic'):
			if(stormObj[ed]['month'][0] in [1, 2, 3, 4, 10, 11, 12]):
				#Calculate the pressure falls over the storm's lifetime
				bergeron = storm.calculate_bergeron(stormObj[ed])
				if(any(b >= 1 for b in bergeron)):
					lon, lat = stormObj[ed]['lon'], stormObj[ed]['lat']				
					lon[lon <  0] = 360 + lon[lon < 0]
					
					plt.plot(lon, lat, 'r-', linewidth=0.5, alpha=0.35, transform=ccrs.PlateCarree())
	# Show the plot.
	plt.title('Storm Tracks')
	plt.savefig('figures/storm_bomb_tracks', bbox_inches='tight', pad_inches=0.05, dpi=300)	
    
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

def plot_timeline(stormObj):
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

	months = [1, 2, 3, 4, 10, 11, 12]
	years = np.arange(1979, 2019, 1)
	cyclGroups = ["Clipper", "Colorado", "GreatBasin", "GulfOfMexico", "EastCoast"]
	
	for year in years:
		for type in cyclGroups:
			meanLat = 0
			meanLon = 0
			count = 0		
			for ed in range(len(stormObj)):
				if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == type and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] == year):
					lon, lat = stormObj[ed]['lon'][0], stormObj[ed]['lat'][0]	
					if(lon < 0):
						lon += 360
					meanLat += lat
					meanLon += lon
					count += 1
			if(count != 0):
				mLat = meanLat / count
				mLon = meanLon / count
				print(type + " (" + str(year) + "): " + str(mLat) + ", " + str(mLon))
				plt.plot(mLon, mLat, linestyle='-', marker='o', color=zoneColors[type], alpha=0.25, markeredgewidth=0, transform=ccrs.PlateCarree())
				plt.text(mLon + 0.3, mLat + 0.3, str(year), alpha=0.5, fontsize=6, transform=ccrs.PlateCarree())
				#ax.annotate('(%s)' % year, xy=(mLon, mLat), textcoords='data', transform=ccrs.PlateCarree())
	plt.title('Cyclogenesis Mean Timeline')
	plt.savefig('figures/cyclogenesis_timeline', bbox_inches='tight', pad_inches=0.05, dpi=300)	

def plot_mean_lat(stormObj):
	fig = plt.figure(figsize=(12,9))

	months = [1, 2, 3, 4, 10, 11, 12]
	years = np.arange(1979, 2019, 1)
	cyclGroups = ["Clipper", "Colorado", "GreatBasin", "GulfOfMexico", "EastCoast"]
	
	for year in years:
		for type in cyclGroups:
			meanLat = 0
			count = 0		
			for ed in range(len(stormObj)):
				if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == type and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] == year):
					lat = stormObj[ed]['lat'][0]	
					meanLat += lat
					count += 1
			if(count != 0):
				mLat = meanLat / count
				plt.plot(year, mLat, '-o', color=zoneColors[type])
	# Add "Fake" Lines for Legend
	for type in cyclGroups:	
		plt.plot(0, 0, linestyle='-', alpha=1, color=zoneColors[type], label=typeStr[type])
	plt.xlim(1979, 2018)
	plt.xlabel("Year")
	plt.ylim(20, 70)
	plt.ylabel("Latitude")
	plt.legend()
	plt.title('Cyclogenesis Mean Latitude')
	plt.savefig('figures/mean_latitude', bbox_inches='tight', pad_inches=0.05, dpi=300)		
	
# List of function calls, uncomment out the plots you want
#plot_positions(stormPosn)
#plot_tracks(storms)
#plot_timeline(storms)
#plot_zones()
#plot_mean_lat(storms)
plot_bomb_tracks(storms)

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