# 07/07/2022
# Vikki Thompson
#
# WASH climate storylines for Nepal


## Add libraries
import os
import cdsapi
import subprocess
import iris
import iris.analysis
import iris.coord_categorisation as icc
from iris.coord_categorisation import add_season_membership
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
import cartopy.feature as cf
plt.ion()
plt.show()

## Define regions
reg = {'Nepal': [23, 33, 78, 91], 'Nepal_west':[23, 33, 78, 85],'Nepal_east':[23, 33, 85, 91], 'Ethiopia':[3, 17, 30, 50], 'Ethiopia_lo':[3, 17, 40, 50], 'Ethiopia_hi':[3, 17, 30, 40], 'Uganda':[-4, 6, 28, 35], 'Uganda_north':[1, 6, 28, 35], 'Uganda_south':[-4, 1, 28, 35], 'Senegal':[11, 18, 341, 349]} # S, N, W, E



def get_region_ts(region, cube_list):
    # for each in cube list extract region
    ts_list = []
    for each in cube_list:
        new_cube = extract_box(each, reg[region])
        if len(np.shape(new_cube.data))== 1:
            pass 
        else:
            ts = new_cube.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)
            ts_list.append(ts)
    return ts_list


def extract_box(cube_in, bdry):
    ' bdry = list, [north, south, east, west] '
    const_lat = iris.Constraint(latitude = lambda cell:bdry[0] < cell < bdry[1])
    const_lon = iris.Constraint(longitude = lambda cell:bdry[2] < cell < bdry[3])
    return cube_in.extract(const_lat & const_lon)


def clear_tmp():
    files = glob.glob('/user/work/hh21501/wash/tmp/*.nc')
    for f in files:
        os.remove(f)
    return

def plot_spagetti(ts_list, ax):
    # plot relative to first decade, and as 10-year smoothed'''
    dec_data = []; dec_year = []
    n=0
    for each in ts_list:
        n+=1
        print(n)
        if len(each.data)<15:
            pass
        else:
            try: 
                icc.add_year(each, each.coord('time'))
            except ValueError:
                pass
            # crop to 2015-2100
            try:
                each = each.extract(iris.Constraint(year=lambda cell: 2014 < cell < 2101))
                each = each.aggregated_by(['year'], iris.analysis.MEAN)
                each.data = each.data - np.mean(each.data[:10])
                # 10-year rolling
                dec = []
                for i in np.arange((len(each.data)-10)):
                    dec.append(np.mean(each.data[i:i+10]))
                # +5 to make centre of decade plotted
                ax.plot(each.coord('year').points[:-10]+5, dec, color='k', alpha=0.2)
                dec_data.append(dec)
                dec_year.append(each.coord('year').points[:-10]+5)
            except AttributeError:
                pass
    # plot ensemble mean
    for n, each in enumerate(dec_data):
        for i in range(dec_year[n][0]-2020):
            each.insert(0, np.nan)
        for i in range(2096-dec_year[n][-1]): # append nan this many times
            each.append(np.nan)  
    for n, each in enumerate(dec_data):
        if len(each) != 77:
            dec_data.remove(each)
    ax.plot(np.arange(2020, 2097, 1), np.nanmean(dec_data, axis = 0), color='r', LineWidth=2)  
    # format plot
    ax.set_xlim([2020, 2095])
    ax.axvspan(2045, 2055, color='red', alpha=0.3)
    ax.axvspan(2085, 2095, color='red', alpha=0.3)
    ax.set_xlabel('Year')
    ax.grid()
    return 


def timeslice(cube_list, y1, y2):
    new_list= iris.cube.CubeList([])
    n=0
    for each in cube_list:
        n+=1
        print(n)
        try: 
            icc.add_year(each, each.coord('time'))
        except ValueError:
            pass
            # crop to 2015-2100
        try:
            each = each.extract(iris.Constraint(year=lambda cell: y1 < cell < y2))
            if len(each.coord('time').points) < (y2-y1-1):
                pass
            else:
                each = each.collapsed('time', iris.analysis.MEAN)
                #each.add_dim_coord('Realization', n)
                new_list.append(each)
        except AttributeError:
            pass
    print('area mean calculation')
    new_cube = new_list[0].copy()
    for i in range(np.shape(new_cube)[0]):
        for j in range(np.shape(new_cube)[1]):
            data_list=[]
            for each in new_list:
                data_list.append(each.data[i,j])
            new_cube.data[i,j]=np.nanmean(data_list)
    return new_cube


def regrid_cubelist(cube_list):
    # identifyhighest res cube
    size_res = []
    for each in cube_list:
        size_res.append(np.shape(each)[1]*np.shape(each)[2])
    new_grid = cube_list[size_res.index(np.max(size_res))]
    new_list = []
    for each in cube_list:
        mod_cs = each.coord_system(iris.coord_systems.CoordSystem)
        new_grid.coord(axis='x').coord_system = mod_cs
        new_grid.coord(axis='y').coord_system = mod_cs
        new_cube = each.regrid(new_grid, iris.analysis.Linear())
        new_list.append(new_cube)
    return new_list


def future_maps(cube_mid, cube_end, crange, region, axs):
    lats=cube_mid.coord('latitude').points
    lons=cube_mid.coord('longitude').points
    c = axs[0].contourf(lons,lats,cube_mid.data, levels=crange, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axs[0].set_title('2041-2060')
    c = axs[1].contourf(lons,lats,cube_end.data, levels=crange, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axs[1].set_title('2081-2100')
    for each in axs:
        each.add_feature(cf.BORDERS)
        each.add_feature(cf.COASTLINE)
        each.set_xlim([region[2]-1, region[3]+1])
        each.set_ylim([region[0]-1, region[1]+2])
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.3, 0.02, 0.4])
    fig.colorbar(c, cax=cbar_ax, ticks=[crange[0], 0, crange[-1]])




### TAS DJF
# Timeseries
tas_djf = iris.load('/user/work/hh21501/wash/tas_ssp242_djf.nc') 
fig, axs = plt.subplots(1, 3, figsize=(8, 3))
# change plot size (so no resizing at end)
nep_list = get_region_ts('Uganda', tas_djf)
#nep_list.remove(nep_list[9]) # remove [9] as years don't align (missing 1950s)
plot_spagetti(nep_list, axs[0])
axs[0].set_ylabel('Temp, *C')
axs[0].set_title('Uganda')
north_list = get_region_ts('Uganda_north', tas_djf)
#nep_list.remove(nep_list[9])
plot_spagetti(north_list, axs[1])
axs[1].set_title('North')
south_list = get_region_ts('Uganda_south', tas_djf)
#south_list.remove(nep_list[9])
plot_spagetti(south_list, axs[2])
axs[2].set_title('South')
for ax in axs:
    ax.set_ylim([-.5, 3])

plt.tight_layout()
plt.savefig('Uganda_tas_djf_ts.png')

# Maps
x_regrid = regrid_cubelist(tas_djf) # all to same grid
x_now = timeslice(x_regrid, 2020, 2041)
x_mid = timeslice(x_regrid, 2041, 2061) # period 1
cube_mid = x_mid - x_now
x_end = timeslice(x_regrid, 2081, 2101) # period 2
cube_end = x_end - x_now
fig, axs = plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,3))
future_maps(cube_mid, cube_end, np.linspace(-2, 2, 20), reg['Uganda'], axs)
for ax in axs:
    ax.plot([29, 35], [1, 1], 'k')

plt.savefig('Uganda_tas_djf_map.png')


### TAS JJA
# Timeseries
tas_djf = iris.load('/user/work/hh21501/wash/tas_ssp242_jja.nc')
fig, axs = plt.subplots(1, 3, figsize=(8, 3))
# change plot size (so no resizing at end)
nep_list = get_region_ts('Uganda', tas_djf)
plot_spagetti(nep_list, axs[0])
axs[0].set_ylabel('Temp, *C')
axs[0].set_title('Uganda')
north_list = get_region_ts('Uganda_north', tas_djf)
plot_spagetti(north_list, axs[1])
axs[1].set_title('North')
south_list = get_region_ts('Uganda_south', tas_djf)
plot_spagetti(south_list, axs[2])
axs[2].set_title('South')
for ax in axs:
    ax.set_ylim([-.5, 3])

plt.tight_layout()
plt.savefig('Uganda_tas_jja_ts.png')

# Maps
x_regrid = regrid_cubelist(tas_djf) # all to same grid
x_now = timeslice(x_regrid, 2020, 2041)
x_mid = timeslice(x_regrid, 2041, 2061) # period 1
cube_mid = x_mid - x_now
x_end = timeslice(x_regrid, 2081, 2101) # period 2
cube_end = x_end - x_now
fig, axs = plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,3))
future_maps(cube_mid, cube_end, np.linspace(-2, 2, 20), reg['Uganda'], axs)
for ax in axs:
    ax.plot([29, 35], [1, 1], 'k')

plt.savefig('Uganda_tas_jja_map.png')



### TASMAX JJA
# Timeseries
tas_djf = iris.load('/user/work/hh21501/wash/tasmax_ssp242_jja.nc')
fig, axs = plt.subplots(1, 3, figsize=(8, 3))
# change plot size (so no resizing at end)
nep_list = get_region_ts('Uganda', tas_djf)
plot_spagetti(nep_list, axs[0])
axs[0].set_ylabel('Temp, *C')
axs[0].set_title('Uganda')
north_list = get_region_ts('Uganda_north', tas_djf)
plot_spagetti(north_list, axs[1])
axs[1].set_title('North')
south_list = get_region_ts('Uganda_south', tas_djf)
plot_spagetti(south_list, axs[2])
axs[2].set_title('South')
for ax in axs:
    ax.set_ylim([-1, 4])

plt.tight_layout()
plt.savefig('Uganda_tasmax_jja_ts.png')

# Maps
x_regrid = regrid_cubelist(tas_djf) # all to same grid
x_now = timeslice(x_regrid, 2020, 2041)
x_mid = timeslice(x_regrid, 2041, 2061) # period 1
cube_mid = x_mid - x_now
x_end = timeslice(x_regrid, 2081, 2101) # period 2
cube_end = x_end - x_now
fig, axs = plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,3))
future_maps(cube_mid, cube_end, np.linspace(-2, 2, 20), reg['Uganda'], axs)
for ax in axs:
    ax.plot([29, 35], [1, 1], 'k')

plt.savefig('Uganda_tasmax_jja_map.png')




### PR JJA
# Timeseries
tas_djf = iris.load('/user/work/hh21501/wash/pr_ssp242_jja.nc')
pr = []
for each in tas_djf:
    pr.append(each*86400*90) # to mm/season (90 days)

fig, axs = plt.subplots(1, 3, figsize=(8, 3))
# change plot size (so no resizing at end)
nep_list = get_region_ts('Uganda', pr)
plot_spagetti(nep_list, axs[0])
axs[0].set_ylabel('Pr, mm')
axs[0].set_title('Uganda')
nep_list = get_region_ts('Uganda_north', pr)
plot_spagetti(nep_list, axs[1])
axs[1].set_title('North')
nep_list = get_region_ts('Uganda_south', pr)
plot_spagetti(nep_list, axs[2])
axs[2].set_title('South')
for ax in axs:
    ax.set_ylim([-100, 200])

plt.tight_layout()
plt.savefig('Uganda_pr_jja_ts.png')

# Maps
x_regrid = regrid_cubelist(pr) # all to same grid
x_now = timeslice(x_regrid, 2020, 2041)
x_mid = timeslice(x_regrid, 2041, 2061) # period 1
cube_mid = x_mid - x_now
x_end = timeslice(x_regrid, 2081, 2101) # period 2
cube_end = x_end - x_now
fig, axs = plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,3))
future_maps(cube_mid, cube_end, np.linspace(-50, 50, 20), reg['Ethiopia'], axs)
for ax in axs:
    ax.plot([29, 35], [1, 1], 'k')

plt.savefig('Uganda_pr_jja_map.png')


### PR MAM
# Timeseries
tas_djf = iris.load('/user/work/hh21501/wash/pr_ssp242_mam.nc')
pr = []
for each in tas_djf:
    pr.append(each*86400) # to mm/season (90 days)

fig, axs = plt.subplots(1, 3, figsize=(8, 3))
# change plot size (so no resizing at end)
nep_list = get_region_ts('Uganda', pr)
plot_spagetti(nep_list, axs[0])
axs[0].set_ylabel('Pr, mm')
axs[0].set_title('Uganda')
north_list = get_region_ts('Uganda_north', pr)
plot_spagetti(north_list, axs[1])
axs[1].set_title('North')
south_list = get_region_ts('Uganda_south', pr)
plot_spagetti(south_list, axs[2])
axs[2].set_title('South')
for ax in axs:
    ax.set_ylim([-50, 100])

plt.tight_layout()
plt.savefig('Uganda_pr_mam_ts.png')

# Maps
x_regrid = regrid_cubelist(pr) # all to same grid
x_now = timeslice(x_regrid, 2020, 2041)
x_mid = timeslice(x_regrid, 2041, 2061) # period 1
cube_mid = x_mid - x_now
x_end = timeslice(x_regrid, 2081, 2101) # period 2
cube_end = x_end - x_now
fig, axs = plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,3))
future_maps(cube_mid, cube_end, np.linspace(-20, 20, 20), reg['Uganda'], axs)
for ax in axs:
    ax.plot([29, 35], [1, 1], 'k')

plt.savefig('Uganda_pr_mam_map.png')


### PR OND
# Timeseries
tas_djf = iris.load('/user/work/hh21501/wash/pr_ssp242_ond.nc')
pr = []
for each in tas_djf:
    pr.append(each*86400) # to mm/season (90 days)

fig, axs = plt.subplots(1, 3, figsize=(8, 3))
# change plot size (so no resizing at end)
nep_list = get_region_ts('Uganda', pr)
plot_spagetti(nep_list, axs[0])
axs[0].set_ylabel('Pr, mm')
axs[0].set_title('Uganda')
north_list = get_region_ts('Uganda_north', pr)
plot_spagetti(north_list, axs[1])
axs[1].set_title('North')
south_list = get_region_ts('Uganda_south', pr)
plot_spagetti(south_list, axs[2])
axs[2].set_title('South')
for ax in axs:
    ax.set_ylim([-50, 100])

plt.tight_layout()
plt.savefig('Uganda_pr_ond_ts.png')

# Maps
x_regrid = regrid_cubelist(pr) # all to same grid
x_now = timeslice(x_regrid, 2020, 2041)
x_mid = timeslice(x_regrid, 2041, 2061) # period 1
cube_mid = x_mid - x_now
x_end = timeslice(x_regrid, 2081, 2101) # period 2
cube_end = x_end - x_now
fig, axs = plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,3))
future_maps(cube_mid, cube_end, np.linspace(-40, 40, 20), reg['Uganda'], axs)
for ax in axs:
    ax.plot([29, 35], [1, 1], 'k')

plt.savefig('Uganda_pr_ond_map.png')



### PRMAX JJA
# Timeseries
tas_djf = iris.load('/user/work/hh21501/wash/pr2_ssp242_jja.nc')
pr = []
for each in tas_djf:
    pr.append(each*86400) # to mm/day

pr = pr[1:]

fig, axs = plt.subplots(1, 3, figsize=(8, 3))
# change plot size (so no resizing at end)
nep_list = get_region_ts('Uganda', pr)
plot_spagetti(nep_list, axs[0])
axs[0].set_ylabel('Pr, mm')
axs[0].set_title('Uganda')
nep_list = get_region_ts('Uganda_north', pr)
plot_spagetti(nep_list, axs[1])
axs[1].set_title('North')
nep_list = get_region_ts('Uganda_south', pr)
plot_spagetti(nep_list, axs[2])
axs[2].set_title('South')
for ax in axs:
    ax.set_ylim([-5, 10])

plt.tight_layout()
plt.savefig('Uganda_pr2_jja_ts.png')

# Maps
x_regrid = regrid_cubelist(pr) # all to same grid
x_now = timeslice(x_regrid, 2020, 2041)
x_mid = timeslice(x_regrid, 2041, 2061) # period 1
cube_mid = x_mid - x_now
x_end = timeslice(x_regrid, 2081, 2101) # period 2
cube_end = x_end - x_now
fig, axs = plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,3))
future_maps(cube_mid, cube_end, np.linspace(-4, 4, 20), reg['Uganda'], axs)
for ax in axs:
    ax.plot([29, 35], [1, 1], 'k')

plt.savefig('Uganda_pr2_jja_map.png')




