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
import scipy
from scipy import spatial
from scipy import stats
import storm_functions as storm

# Load storm data
det_storms = np.load('storm_det_slp.npz', encoding='latin1', allow_pickle=True)
stormPosn = det_storms['storms']

tracked_storms = np.load('storm_track_slp.npz', encoding='latin1', allow_pickle=True)
storms = tracked_storms['storms']

census_storms = np.load('storm_census.npz', encoding='latin1', allow_pickle=True)

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

def get_projection_object(type="Lambert"):
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
	if type == "Lambert":
		projObj = ccrs.LambertConformal(central_latitude=cen_lat, 
										central_longitude=cen_lon,
										standard_parallels=std_pll,
										globe=globe,
										cutoff=cutoff)
	elif type == "EqualArea":
		projObj = ccrs.AlbersEqualArea(central_latitude=cen_lat,
									   central_longitude=cen_lon,
									   standard_parallels=std_pll,
									   globe=globe)
	else:
		print("get_projection_object(): Invalid type " + type + " sent to function")
		return None
	return projObj

def plot_positions(stormObj):
	# Grab the projection object from the netCDF file
	projObj = get_projection_object()    
	# Create our figure and axis object
	fig = plt.figure(figsize=(12,9))
	ax = plt.axes(projection=projObj)

	ax.set_extent(plotExtent, crs=ccrs.PlateCarree())   

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
	
def plot_zones(stormObj = None):
	print("plot_zones()")
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
	
	if stormObj is not None:
		print("plot_zones(): stormObj is not None")
		for ed in range(len(stormObj)):
			if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == 'Other'):	
				if (stormObj[ed]['month'][0] in [1, 2, 3, 4, 10, 11, 12]):
					lon, lat = stormObj[ed]['lon'], stormObj[ed]['lat']
					lon[lon <  0] = 360 + lon[lon < 0]
					plt.plot(lon[0], lat[0], 'bo', alpha=0.25, markeredgewidth=0, markersize=2, transform=ccrs.PlateCarree())
					
		plt.legend(framealpha = 1.0)
		plt.title('Storm Classification Zones (Other Marked)')
		plt.savefig('figures/storm_classes_with_other', bbox_inches='tight', pad_inches=0.05, dpi=300)		
	else:
		print("plot_zones(): stormObj is None")
		plt.legend(framealpha = 1.0)
		plt.title('Storm Classification Zones')
		plt.savefig('figures/storm_classes', bbox_inches='tight', pad_inches=0.05, dpi=300)	
	print("plot_zones(): Done")
	
def plot_tracks(stormObj):
	print("plot_tracks()")
	projObj = get_projection_object()    
	# Create our figure and axis object
	fig = plt.figure(figsize=(12,9))
	ax = plt.axes(projection=projObj)
	#ax = plt.axes(projection=ccrs.PlateCarree())
	ax.set_extent(plotExtent, crs=ccrs.PlateCarree())   

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
	plt.legend(framealpha = 1.0)
	plt.title('Storm Tracks')
	plt.savefig('figures/storm_tracks', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_tracks(): Done")
	
def plot_bomb_tracks(stormObj):
	print("plot_bomb_tracks()")
	projObj = get_projection_object()    
	# Create our figure and axis object
	fig = plt.figure(figsize=(12,9))
	ax = plt.axes(projection=projObj)
	#ax = plt.axes(projection=ccrs.PlateCarree())
	ax.set_extent(plotExtent, crs=ccrs.PlateCarree())   

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
	print("plot_bomb_tracks(): Done")
	
def plot_mean_track(stormObj, bombOnly = False):
	print("plot_mean_track(bombOnly = " + str(bombOnly) + ")")
	projObj = get_projection_object()    
	# Create our figure and axis object
	fig = plt.figure(figsize=(12,9))
	ax = plt.axes(projection=projObj)
	#ax = plt.axes(projection=ccrs.PlateCarree())
	ax.set_extent(plotExtent, crs=ccrs.PlateCarree())   

	# Draw our plot, coastlines first, then the contours.
	states = cartopy.feature.NaturalEarthFeature(category='cultural',
												name='admin_1_states_provinces_lakes',
												scale='110m',
												facecolor='none')
	ax.add_feature(states, edgecolor='k')    
	ax.coastlines()
	
	# Process Variables
	monthList = [1, 2, 3, 4, 10, 11, 12]
	yearList = np.arange(1979, 2020, 1)	
	cycloneList = ['Clipper', 'Northwest', 'Colorado', 'GreatBasin', 'GulfOfMexico', 'EastCoast']	
		
	for classification in cycloneList:	
		# Find the longest cyclone.
		longest = 0
		for ed in range(len(stormObj)):
			if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['month'][0] in monthList
				and stormObj[ed]['year'][0] in yearList and stormObj[ed]['classification'] == classification):
				count = len(stormObj[ed]['amp'])
				if(count > longest):
					longest = count		
		# Test classification
		meanLat = np.zeros((longest))
		meanLon = np.zeros((longest))
		count = np.zeros((longest))	
		for ed in range(len(stormObj)):
			if (stormObj[ed]['type'] == 'cyclonic'):
				# Test months
				if(storms[ed]['month'][0] in monthList):
					# Test years
					if(storms[ed]['year'][0] in yearList):										
						if(storms[ed]['classification'] == classification):
							bergeron = storm.calculate_bergeron(stormObj[ed])
							if(bombOnly == True):
								if(any(b >= 1 for b in bergeron)):	
									for i in range(len(stormObj[ed]['lat'])):
										meanLat[i] += stormObj[ed]['lat'][i] 
										meanLon[i] += stormObj[ed]['lon'][i]
										count[i] += 1
							else:
								for i in range(len(stormObj[ed]['lat'])):
									meanLat[i] += stormObj[ed]['lat'][i] 
									meanLon[i] += stormObj[ed]['lon'][i]
									count[i] += 1
		meanLat[:] = meanLat[:] / count[:]
		meanLon[:] = meanLon[:] / count[:]
		
		print(classification + ": stdLat: %.3f" % np.std(meanLat))
		
		#print(meanLat)
		#print(meanLon)
		#print(count)
		# Draw the lines.
		plt.plot(meanLon, meanLat, linestyle='-', color=zoneColors[classification], alpha=1, markeredgewidth=0, transform=ccrs.PlateCarree())
	
	legendVals = list(typeStr.values())
	legendVals.remove("Other Lows")
	
	plt.legend(legendVals, loc=0, framealpha = 1.0)
	if(bombOnly == True):
		plt.title("Mean Bomb Cyclone Track (" + str(yearList[0]) + " - " + str(yearList[-1]) + ")")
		plt.savefig('figures/mean_track_bomb', bbox_inches='tight', pad_inches=0.05, dpi=300)
	else:
		plt.title("Mean Cyclone Track (" + str(yearList[0]) + " - " + str(yearList[-1]) + ")")
		plt.savefig('figures/mean_track', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_mean_track(): Done.")
	
def plot_bomb_frequency(stormObj):
	print("plot_bomb_frequency()")
	fig = plt.figure(figsize=(12,9))
	yrs = np.arange(1979, 2020, 1)
	counts = []
	for y in yrs:
		indCount = 0
		for ed in range(len(stormObj)):
			if (stormObj[ed]['type'] == 'cyclonic'):
				if(stormObj[ed]['month'][0] in [1, 2, 3, 4, 10, 11, 12] and stormObj[ed]['year'][0] == y):
					#Calculate the pressure falls over the storm's lifetime
					bergeron = storm.calculate_bergeron(stormObj[ed])
					if(any(b >= 1 for b in bergeron)):
						indCount += 1
		counts.append(indCount)
	#Draw the plot
	plt.bar(yrs, counts)
	plt.xlabel("Year")
	plt.ylabel("Number of Cyclones")
	plt.title("Number of Bomb Cyclones")
	plt.savefig('figures/bomb_frequency', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_bomb_frequency(): Done.")
    
def plot_frequencies(stormCensus):
	print("plot_frequencies()")

	projObj = get_projection_object(type="EqualArea")
	plotExtent = [-140, -60, 20, 70]

	lonCount = ((plotExtent[1]+360) - (plotExtent[0]+360))
	latCount = ((plotExtent[3]) - (plotExtent[2]))
	lon = np.linspace(plotExtent[0], plotExtent[1], num=lonCount)
	lat = np.linspace(plotExtent[2], plotExtent[3], num=latCount)	
	
	llon, llat = np.meshgrid(lon, lat)
	
	fig, axs = plt.subplots(1, 2)
	ax1 = axs[0]
	ax2 = axs[1]
	
	# Distribution of cyclone tracks and intensity (and for anticyclones)
	print("plot_frequencies(): storm_track_distribution")
	#fig = plt.figure(figsize=(12,9))
	ax1 = plt.axes(projection=projObj)
	ax1.set_extent(plotExtent, crs=ccrs.PlateCarree())   
	states = cartopy.feature.NaturalEarthFeature(category='cultural',
												name='admin_1_states_provinces_lakes',
												scale='110m',
												facecolor='none')
	ax1.add_feature(states, edgecolor='k')    
	ax1.coastlines()
	PCM1 = ax1.pcolormesh(llon, llat, stormCensus['count'], transform=ccrs.PlateCarree())
	fig.colorbar(PCM1, ax=ax1)
	
	#H = plt.colorbar()
	#H.set_label('Count')
	#plt.clim(0, 250)
	#plt.savefig('figures/storm_track_distribution.png', bbox_inches='tight', pad_inches=0.05, dpi=300)

	print("plot_frequencies(): genesis_distribution")
	#fig = plt.figure(figsize=(12,9))
	ax2 = plt.axes(projection=projObj)
	ax2.set_extent(plotExtent, crs=ccrs.PlateCarree())   
	ax2.add_feature(states, edgecolor='k')    
	ax2.coastlines()
	PCM2 = ax2.pcolormesh(llon, llat, stormCensus['genesis'], transform=ccrs.PlateCarree())
	fig.colorbar(PCM2, ax=ax2)
	#plt.title('Cyclogenesis Locations')
	#H = plt.colorbar()
	#H.set_label('Count')
	#plt.clim(0, 250)
	plt.savefig('figures/panel_test.png', bbox_inches='tight', pad_inches=0.05, dpi=300)

	"""
	print("plot_frequencies(): termination_distribution")
	fig = plt.figure(figsize=(12,9))
	ax = plt.axes(projection=projObj)
	ax.set_extent(plotExtent, crs=ccrs.PlateCarree())   
	states = cartopy.feature.NaturalEarthFeature(category='cultural',
												name='admin_1_states_provinces_lakes',
												scale='110m',
												facecolor='none')
	ax.add_feature(states, edgecolor='k')    
	ax.coastlines()
	plt.pcolormesh(llon, llat, stormCensus['termination'], transform=ccrs.PlateCarree())
	plt.title('Cyclolysis Locations')
	H = plt.colorbar()
	H.set_label('Count')
	plt.clim(0, 75)
	plt.savefig('figures/termination_distribution.png', bbox_inches='tight', pad_inches=0.05, dpi=300)
	"""
	print("plot_frequencies(): Done")

def plot_timeline(stormObj):
	print("plot_timeline()")
	# Grab the projection object from the netCDF file
	projObj = get_projection_object()    
	# Create our figure and axis object
	fig = plt.figure(figsize=(12,9))
	ax = plt.axes(projection=projObj)
	#ax = plt.axes(projection=ccrs.PlateCarree())
	ax.set_extent(plotExtent, crs=ccrs.PlateCarree())   

	# Draw our plot, coastlines first, then the contours.
	states = cartopy.feature.NaturalEarthFeature(category='cultural',
												name='admin_1_states_provinces_lakes',
												scale='110m',
												facecolor='none')
	ax.add_feature(states, edgecolor='k')    
	ax.coastlines()

	months = [1, 2, 3, 4, 10, 11, 12]
	years = np.arange(1979, 2020, 1)
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
				
				plt.plot(mLon, mLat, linestyle='-', marker='o', color=zoneColors[type], alpha=0.25, markeredgewidth=0, transform=ccrs.PlateCarree())
				plt.text(mLon + 0.3, mLat + 0.3, str(year), alpha=0.5, fontsize=6, transform=ccrs.PlateCarree())
				#ax.annotate('(%s)' % year, xy=(mLon, mLat), textcoords='data', transform=ccrs.PlateCarree())
	plt.title('Cyclogenesis Mean Timeline')
	plt.savefig('figures/cyclogenesis_timeline', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_timeline(): Done")
	
def plot_bomb_frequency_typed(stormObj):
	print("plot_bomb_frequency_typed()")
	fig = plt.figure(figsize=(12,9))
	yrs = np.arange(1979, 2020, 1)
	totals_year = np.zeros((len(yrs)))
	for className in typeStr.keys():
		count_list_ind = np.zeros((len(yrs)))
		for y in yrs:
			indCount = 0
			for ed in range(len(stormObj)):
				if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == className):
					if(stormObj[ed]['month'][0] in [1, 2, 3, 4, 10, 11, 12] and stormObj[ed]['year'][0] == y):
						#Calculate the pressure falls over the storm's lifetime
						bergeron = storm.calculate_bergeron(stormObj[ed])
						if(any(b >= 1 for b in bergeron)):
							indCount += 1
			count_list_ind[y - 1979] = indCount
		#Draw the plot
		plt.bar(yrs, count_list_ind, bottom=totals_year, color=zoneColors[className])
		totals_year += count_list_ind
	plt.ylim(0, 60)
	plt.xlabel("Year")
	plt.ylabel("Number of Cyclones")	
	plt.legend(typeStr.values(), loc=0, framealpha = 1.0)
	plt.title("Number of Bomb Cyclones by Cyclone Type")
	plt.savefig('figures/bomb_frequency_typed', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_bomb_frequency_typed(): Done")
	
def plot_frequency_typed(stormObj):
	print("plot_frequency_typed()")
	fig = plt.figure(figsize=(12,9))
	yrs = np.arange(1979, 2020, 1)
	totals_year = np.zeros((len(yrs)))
	for className in typeStr.keys():
		count_list_ind = np.zeros((len(yrs)))
		for y in yrs:
			indCount = 0
			for ed in range(len(stormObj)):
				if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == className and className != "Other"):
					if(stormObj[ed]['month'][0] in [1, 2, 3, 4, 10, 11, 12] and stormObj[ed]['year'][0] == y):
						indCount += 1
			count_list_ind[y - 1979] = indCount
		#Draw the plot
		plt.bar(yrs, count_list_ind, bottom=totals_year, color=zoneColors[className])
		totals_year += count_list_ind
	
	legendList = list(typeStr.values())
	legendList.remove("Other Lows")
	
	plt.legend(legendList, loc=0, framealpha = 1.0)
	plt.xlabel("Year")
	plt.ylabel("Number of Cyclones")	
	plt.title("Number of Cyclones by Cyclone Type")
	plt.savefig('figures/frequency_typed', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_frequency_typed(): Done")
	
def plot_frequency_for_data(stormObj):
	print("plot_frequency_for_data()")
	fig = plt.figure(figsize=(12,9))
	yrs = np.arange(1979, 2020, 1)
	totals_year = np.zeros((len(yrs)))
	okList = []
	for className in typeStr.keys():
		if className == 'Other':
			continue
		count_list_ind = np.zeros((len(yrs)))
		for y in yrs:
			indCount = 0
			for ed in range(len(stormObj)):
				if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == className):
					if(stormObj[ed]['month'][0] in [1, 2, 3, 4, 10, 11, 12] and stormObj[ed]['year'][0] == y):
						indCount += 1
			count_list_ind[y - 1979] = indCount
		if not (0 in count_list_ind):
			okList.append(className)
	# Redo counts with only non-zeros
	for className in okList:
		count_list_ind = np.zeros((len(yrs)))
		for y in yrs:
			indCount = 0
			for ed in range(len(stormObj)):
				if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == className):
					if(stormObj[ed]['month'][0] in [1, 2, 3, 4, 10, 11, 12] and stormObj[ed]['year'][0] == y):
						indCount += 1
			count_list_ind[y - 1979] = indCount
		plt.bar(yrs, count_list_ind, bottom=totals_year, color=zoneColors[className])
		totals_year += count_list_ind
	plt.legend(typeStr.values(), loc=0, framealpha = 1.0)
	plt.xlabel("Year")
	plt.ylabel("Number of Cyclones")	
	plt.title("Number of Cyclones by Cyclone Type")
	plt.savefig('figures/frequency_typed_non_zeros', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_frequency_for_data(): Done")
	
def plot_strength_bars(stormObj):
	print("plot_strength_bars()")
	fig = plt.figure(figsize=(12,9))

	months = [1, 2, 3, 4, 10, 11, 12]
	years = np.arange(1979, 2020, 1)
	dt = 3
	
	# Find the longest cyclone.
	longest = 0
	for ed in range(len(stormObj)):
		if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] in years):
			count = len(stormObj[ed]['amp'])
			if(count > longest):
				longest = count
				
	totals_column = np.zeros((longest))
	for className in typeStr.keys():
		fall_array = np.zeros((longest))
		for ed in range(len(stormObj)):
			if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] in years and stormObj[ed]['classification'] == className):
				dPressure = np.diff(stormObj[ed]['amp'])
				largestFall = np.nanmin(dPressure)
				index_of_largest_fall = np.argmin(dPressure)
				time_of_fall = index_of_largest_fall * dt
				fall_array[index_of_largest_fall] += 1
		#idx_of_last_non_zero = np.max(np.nonzero(fall_array))
		#splitArray = fall_array[0:idx_of_last_non_zero+1]
		times = np.arange(0, (len(fall_array)*3), 3)
	
		plt.bar(times+1, fall_array, bottom=totals_column, width=2, color=zoneColors[className])
		totals_column += fall_array
	idx_of_last_non_zero = np.max(np.nonzero(fall_array))
	plt.xlim(0, idx_of_last_non_zero*3)
	plt.legend(typeStr.values(), loc=0, framealpha = 1.0)
	plt.xlabel("Time after cyclogenesis (hours)")
	plt.ylabel("Count")
	plt.yscale('log')
	plt.title("Time of Maximum Pressure Fall")
	plt.savefig('figures/pressure_fall_occurance', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_strength_bars(): Done")
	
def plot_strength_bars_bombs(stormObj):
	print("plot_strength_bars_bombs()")
	fig = plt.figure(figsize=(12,9))

	months = [1, 2, 3, 4, 10, 11, 12]
	years = np.arange(1979, 2020, 1)
	dt = 3
	
	# Find the longest cyclone.
	longest = 0
	for ed in range(len(stormObj)):
		if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] in years):
			count = len(stormObj[ed]['amp'])
			if(count > longest):
				longest = count
	
	totals_column = np.zeros((longest))
	for className in typeStr.keys():
		fall_array = np.zeros((longest))
		for ed in range(len(stormObj)):
			if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] in years and stormObj[ed]['classification'] == className):
				bergeron = storm.calculate_bergeron(stormObj[ed])
				if(any(b >= 1 for b in bergeron)):		
					dPressure = np.diff(stormObj[ed]['amp'])
					largestFall = np.nanmin(dPressure)
					index_of_largest_fall = np.argmin(dPressure)
					time_of_fall = index_of_largest_fall * dt
					fall_array[index_of_largest_fall] += 1
		times = np.arange(0, (len(fall_array)*3), 3)
		plt.bar(times+1, fall_array, bottom=totals_column, width=2, color=zoneColors[className])
		totals_column += fall_array
	idx_of_last_non_zero = np.max(np.nonzero(fall_array))
	plt.xlim(0, idx_of_last_non_zero*3)
	plt.legend(typeStr.values(), loc=0, framealpha = 1.0)		
	plt.xlabel("Time after cyclogenesis (hours)")
	plt.ylabel("Count")
	plt.title("Time of Maximum 3-Hour Pressure Fall (Bomb Cyclones Only)")
	plt.savefig('figures/pressure_fall_bombs_occurance', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_strength_bars_bombs(): Done")
	
def plot_cyclone_lifespan(stormObj, bombOnly = False):
	print("plot_cyclone_lifespan(bombOnly = " + str(bombOnly) + ")")
	fig = plt.figure(figsize=(12,9))

	months = [1, 2, 3, 4, 10, 11, 12]
	years = np.arange(1979, 2020, 1)
	dt = 3	
	
	# Find the longest cyclone.
	longest = 0
	for ed in range(len(stormObj)):
		if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] in years):
			count = len(stormObj[ed]['amp'])
			if(count > longest):
				longest = count
	
	lifespan_array = np.zeros((longest))

	for ed in range(len(stormObj)):
		if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] in years):
			bergeron = storm.calculate_bergeron(stormObj[ed])	
			if(bombOnly == True):
				if(any(b >= 1 for b in bergeron)):
					lifespan = len(stormObj[ed]['amp'])
					lifespan_array[lifespan-1] += 1
			else:
				lifespan = len(stormObj[ed]['amp'])
				lifespan_array[lifespan-1] += 1			
				
	idx_of_last_non_zero = np.max(np.nonzero(lifespan_array))
	splitArray = lifespan_array[0:idx_of_last_non_zero+1]
	times = np.arange(0, (len(splitArray)*3), 3)				
				
	plt.bar(times, splitArray)
	plt.xlabel("Cyclone Lifespan (Hours)")
	plt.ylabel("Count")
	if(bombOnly == True):
		plt.title("Lifespan of Cyclones (Bomb Cyclones Only)")
		plt.savefig('figures/lifespan_bombs', bbox_inches='tight', pad_inches=0.05, dpi=300)		
	else:
		plt.yscale('log')
		plt.title("Lifespan of Cyclones")
		plt.savefig('figures/lifespan', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_cyclone_lifespan(): Done")

def plot_mean_lat(stormObj):
	print("plot_mean_lat()")
	fig = plt.figure(figsize=(12,9))

	months = [1, 2, 3, 4, 10, 11, 12]
	years = np.arange(1979, 2020, 1)
	cyclGroups = ["Clipper", "Colorado", "Northwest", "GreatBasin", "GulfOfMexico", "EastCoast"]
	
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
	# Regression Lines
	for type in cyclGroups:
		latList = []
		for year in years:
			lat = 0
			meanLat = 0
			mLat = 0
			count = 0
			for ed in range(len(stormObj)):
				if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == type and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] == year):
					lat = stormObj[ed]['lat'][0]
					meanLat += lat
					count += 1
			if count != 0:
				mLat = meanLat / count
			else:
				mLat = np.nan
			
			latList.append(mLat)
			print(type + " (" + str(year) + "): " + str(mLat))
			
		finiteYmask = np.isfinite(np.array(latList))
		Yclean = np.array(latList)[finiteYmask]
		Xclean = np.array(years)[finiteYmask]			

		slope, intercept, r_value, p_value, std_err = stats.linregress(Xclean, Yclean)
		line = slope*years+intercept	
		plt.plot(years, line, color=zoneColors[type], linestyle='--')
		
		print("Statistics " + type + ": m: %.3f, b: %.3f, r: %.3f, p: %.3f, stdErr: %.3f, stdDev: %.3f" % (slope, intercept, r_value, p_value, std_err, np.std(Yclean)))
	# Add "Fake" Lines for Legend
	for type in cyclGroups:	
		plt.plot(0, 0, linestyle='-', alpha=1, color=zoneColors[type], label=typeStr[type])
	plt.xlim(1979, 2020)
	plt.xlabel("Year")
	plt.ylim(20, 70)
	plt.ylabel("Latitude")
	plt.legend(framealpha = 1.0)
	plt.title('Cyclogenesis Mean Latitude')
	plt.savefig('figures/mean_latitude', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_mean_lat(): Done")	
	
def plot_mean_pressure(stormObj):
	print("plot_mean_pressure()")
	fig = plt.figure(figsize=(12,9))

	months = [1, 2, 3, 4, 10, 11, 12]
	years = np.arange(1979, 2020, 1)
	cyclGroups = ["Clipper", "Colorado", "Northwest", "GreatBasin", "GulfOfMexico", "EastCoast"]
	
	for year in years:
		for type in cyclGroups:
			meanPressure = 0
			count = 0		
			for ed in range(len(stormObj)):
				if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == type and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] == year):
					lowestP = np.nanmin(stormObj[ed]['amp']) / 100
					meanPressure += lowestP
					count += 1
			if(count != 0):
				mP = meanPressure / count
				plt.plot(year, mP, '-o', color=zoneColors[type])	
	
	# Regression Lines
	for type in cyclGroups:
		PList = []
		for year in years:
			pres = 0
			meanPres = 0
			mP = 0
			count = 0
			for ed in range(len(stormObj)):
				if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == type and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] == year):
					lowestP = np.nanmin(stormObj[ed]['amp']) / 100
					meanPres += lowestP
					count += 1
			if count != 0:
				mP = meanPres / count
			else:
				mP = np.nan
			
			PList.append(mP)
			print(type + " (" + str(year) + "): " + str(mP))
			
		finiteYmask = np.isfinite(np.array(PList))
		Yclean = np.array(PList)[finiteYmask]
		Xclean = np.array(years)[finiteYmask]			

		slope, intercept, r_value, p_value, std_err = stats.linregress(Xclean, Yclean)
		line = slope*years+intercept	
		plt.plot(years, line, color=zoneColors[type], linestyle='--')
		
		print("Statistics " + type + ": m: %.3f, b: %.3f, r: %.3f, p: %.3f, stdErr: %.3f, stdDev: %.3f" % (slope, intercept, r_value, p_value, std_err, np.std(Yclean)))
	# Add "Fake" Lines for Legend
	for type in cyclGroups:	
		plt.plot(0, 0, linestyle='-', alpha=1, color=zoneColors[type], label=typeStr[type])
	plt.xlim(1979, 2020)
	plt.ylim(950, 1020)
	plt.xlabel("Year")
	plt.ylabel("Pressure (hPa)")
	plt.legend(framealpha = 1.0)
	plt.title('Cyclone Mean Lowest Pressure')
	plt.savefig('figures/mean_lowest_pressure', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_mean_pressure(): Done")

def plot_mean_pressure_bombs(stormObj):
	print("plot_mean_pressure_bombs()")
	fig = plt.figure(figsize=(12,9))

	months = [1, 2, 3, 4, 10, 11, 12]
	years = np.arange(1979, 2020, 1)
	cyclGroups = ["Clipper", "Colorado", "GreatBasin", "GulfOfMexico", "EastCoast"]
	
	for year in years:
		for type in cyclGroups:
			meanPressure = 0
			count = 0		
			for ed in range(len(stormObj)):
				if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['classification'] == type and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] == year):
					bergeron = storm.calculate_bergeron(stormObj[ed])
					if(any(b >= 1 for b in bergeron)):
						#lat = stormObj[ed]['lat'][0]	
						lowestP = np.nanmin(stormObj[ed]['amp'])
						meanPressure += lowestP
						count += 1
			if(count != 0):
				meanPressure = (meanPressure / count) / 100
				plt.plot(year, meanPressure, '-o', color=zoneColors[type])			
	# Add "Fake" Lines for Legend
	for type in cyclGroups:	
		plt.plot(0, 0, linestyle='-', alpha=1, color=zoneColors[type], label=typeStr[type])
	plt.xlim(1979, 2020)
	plt.ylim(900, 1000)
	plt.xlabel("Year")
	plt.ylabel("Pressure")
	plt.legend(framealpha = 1.0)
	plt.title('Cyclone Mean Lowest Pressure (Bomb Cyclones Only)')
	plt.savefig('figures/mean_lowest_pressure_bomb', bbox_inches='tight', pad_inches=0.05, dpi=300)
	print("plot_mean_pressure_bombs(): Done")
	
def plot_lifespan_panel(stormObj):
	print("plot_lifespan_panel()")
	
	fig, axs = plt.subplots(1, 2)
	ax1 = axs[0]
	ax2 = axs[1]
	
	months = [1, 2, 3, 4, 10, 11, 12]
	years = np.arange(1979, 2020, 1)
	dt = 3	
	
	longest = 0
	for ed in range(len(stormObj)):
		if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] in years):
			count = len(stormObj[ed]['amp'])
			if(count > longest):
				longest = count
	
	# Panel A
	lifespan_array1 = np.zeros((longest))
	for ed in range(len(stormObj)):
		if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] in years):
			lifespan = len(stormObj[ed]['amp'])
			lifespan_array1[lifespan-1] += 1			
				
	idx_of_last_non_zero = np.max(np.nonzero(lifespan_array1))
	splitArray1 = lifespan_array1[0:idx_of_last_non_zero+1]
	times1 = np.arange(0, (len(splitArray1)*3), 3)				
				
	#plt.bar(times, splitArray)
	ax1.scatter(times1, splitArray1, s=2)
	ax1.set(xlabel='Cyclone Lifespan (Hours)', ylabel='Count')
	
	#ax1.
	
	# Panel B
	lifespan_array2 = np.zeros((longest))
	for ed in range(len(stormObj)):
		if (stormObj[ed]['type'] == 'cyclonic' and stormObj[ed]['month'][0] in months and stormObj[ed]['year'][0] in years):
			bergeron = storm.calculate_bergeron(stormObj[ed])	
			if(any(b >= 1 for b in bergeron)):			
				lifespan = len(stormObj[ed]['amp'])
				lifespan_array2[lifespan-1] += 1			
				
	idx_of_last_non_zero = np.max(np.nonzero(lifespan_array2))
	splitArray2 = lifespan_array2[0:idx_of_last_non_zero+1]
	times2 = np.arange(0, (len(splitArray2)*3), 3)				
				
	ax2.bar(times2, splitArray2)
	#ax1.scatter(times2, splitArray2)
	ax2.set(xlabel='Cyclone Lifespan (Hours)', ylabel='Count')	
	
	plt.suptitle('Lifespan of Cyclones')
	plt.savefig('figures/lifespan_panel', bbox_inches='tight', pad_inches=0.05, dpi=300)	
	print("plot_lifespan_panel(): Done")
	
	
print("Calling plotting routines")
# List of function calls, uncomment out the plots you want
#plot_positions(stormPosn)
#plot_frequencies(census_storms)
#plot_tracks(storms)
#plot_timeline(storms)
#plot_zones()
#plot_mean_lat(storms)
#plot_bomb_tracks(storms)
#plot_bomb_frequency(storms)
#plot_bomb_frequency_typed(storms)
#plot_frequency_typed(storms)
#plot_frequency(storms, months=[1,2,3,4,10,11,12], years=[1980], cyclones_to_include = ['Clipper'])
#plot_strength_bars(storms)
#plot_strength_bars_bombs(storms)
#plot_cyclone_lifespan(storms, bombOnly = False)
#plot_cyclone_lifespan(storms, bombOnly = True)
#plot_mean_track(storms, bombOnly = False)
#plot_mean_track(storms, bombOnly = True)
#plot_mean_pressure(storms)
#plot_mean_pressure_bombs(storms)
#plot_frequency_for_data(storms)
#plot_zones(storms)
plot_lifespan_panel(storms)

print("Done.")