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
data = np.load('storm_track_slp.npz')
storms = data['storms']

# Settings
plotExtent = [-125, -70, 23, 49] # Where to show the plot (LonMin, LonMax, LatMin, LatMax)

def get_projection_object(ncFile):
    # Cartopy has a "globe" object to define more projection standards, we create this first
    globe = ccrs.Globe(ellipse=None,
                       semimajor_axis=6370000,
                       semiminor_axis=6370000,
                       nadgrids="@null")    
    # Now we can create the projection object
    cen_lat  = float(ncFile.centerlat)
    cen_lon  = float(ncFile.centerlon)
    std_pll  = [ncFile.standardpar1]
    cutoff   = -30.0
    projObj = ccrs.LambertConformal(central_latitude=cen_lat, 
                                    central_longitude=cen_lon,
                                    standard_parallels=std_pll,
                                    globe=globe,
                                    cutoff=cutoff)
    
    return projObj

def plot_sample(ncFile):
    # Grab the projection object from the netCDF file
    projObj = get_projection_object(ncFile)    
    # Create our figure and axis object
    fig = plt.figure(figsize=(12,9))
    ax = plt.axes(projection=projObj)
    # Note the coordinate flip here, this is because WRF stores N_S first, then E_W
    lats = ncFile.variables['lat'][:] 
    lons = ncFile.variables['lon'][:] 
       
    clevs = np.arange(980, 1060, 1)
    clevs2 = np.arange(980, 1060, 4)
    
    # Grab the surface temperature from the netCDF file (0th time step, 0th vertical level)
    MSLP = ncFile.variables["prmsl"][0] / 100
    
    ax.set_extent([-125, -70, 23, 49], crs=ccrs.PlateCarree())   
    
    # Draw our plot, coastlines first, then the contours.
    states = cartopy.feature.NaturalEarthFeature(category='cultural',
                                                name='admin_1_states_provinces_lakes',
                                                scale='110m',
                                                facecolor='none')
    ax.add_feature(states, edgecolor='k')    
    ax.coastlines()
    
    plt.contour(lons, lats, MSLP, clevs2, colors='k', linewidths=1, linestyles='dotted', transform=ccrs.PlateCarree()) #, alpha=0.5, transform=ccrs.PlateCarree())
    contours = plt.contourf(lons, lats, MSLP, clevs, cmap="jet", transform=ccrs.PlateCarree())
    plt.colorbar(ax=ax, orientation="horizontal", pad=.05)
    # Show the plot.
    plt.show()
    
plot_sample(nc1)

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
