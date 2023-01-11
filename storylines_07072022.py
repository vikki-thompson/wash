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

## Load data
# Ensemble mean from a selection of CMIP6 (regional) models


def seasonal_timeseries(infile, seas, label):
    ''' seas = 'JJA'  years = 'XXXX' '''
    ''' '''
    # Merge to one file 2015-2100
    cdo_cmd = ['cdo', 'mergetime', infile, 'tmp/ssp'+label+'_p1.nc']
    ret = subprocess.call(cdo_cmd)
    # Extract JJA months
    cdo_seas_p1 = ['cdo', 'selseas,'+seas, 'tmp/ssp'+label+'_p1.nc', 'tmp/ssp'+label+'_p2.nc']
    ret = subprocess.call(cdo_seas_p1)
    # JJA annual values
    cdo_seas_p2 = ['cdo', 'seasmean', 'tmp/ssp'+label+'_p2.nc', 'tmp/ssp_'+label+'.nc']
    ret = subprocess.call(cdo_seas_p2)
    return iris.load('tmp/ssp_'+label+'.nc')[0]


def apply_nepal_mask(cube):
    new_cube = cube.copy()
    NPL=gpd.read_file('nepal_shape/hermes_NPL_new_wgs_0.shp')
    #get the latitude-longitude grid from netcdf file
    lats=new_cube.coord('latitude').points
    lons=new_cube.coord('longitude').points
    lon2d, lat2d = np.meshgrid(lons, lats)
    # reshape 
    lon2 = lon2d.reshape(-1)
    lat2 = lat2d.reshape(-1)
    # 
    mask = []
    for lat, lon in zip(lat2, lon2):
        this_point = gpd.GeoSeries([Point(lon, lat)])
        res = NPL.geometry.contains(this_point)
        mask.append(res.values[0])
    mask = np.array(mask).reshape(lon2d.shape)
    mask = ~mask
    #dim_map = (cube.coord_dims('longitude')[0], cube.coord_dims('latitude')[0])
    #cube_mask = iris.util.broadcast_to_shape(mask, cube.shape, dim_map)
    data = cube.data
    masked_data = np.ma.masked_array(data, mask)
    new_cube.data = masked_data
    return new_cube

def extract_box(cube_in, bdry):
    ' bdry = list, [north, south, east, west] '
    const_lat = iris.Constraint(latitude = lambda cell:bdry[0] < cell < bdry[1])
    const_lon = iris.Constraint(longitude = lambda cell:bdry[2] < cell < bdry[3])
    new_cube = cube_in.extract(const_lat & const_lon)
    return new_cube

def multiple_maps(cubes, levels, cmap):
    fig, ax = plt.subplots(nrows=3, ncols=3, subplot_kw={'projection': ccrs.PlateCarree()})
    ax_list = [ax[0,0], ax[0,1], ax[0,2], ax[1,0], ax[1,1], ax[1,2], ax[2,0], ax[2,1], ax[2,2]]
    for i, each in enumerate(cubes):
        c = map_plot(each, ax_list[i], levels, cmap)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(c, cax=cbar_ax)
    ax[0,0].set_title('2021-2040')
    ax[0,1].set_title('2041-2060')
    ax[0,2].set_title('2081-2100')
    ax[0,0].set_ylabel('SSP1-26')
    ax[1,0].set_ylabel('SSP2-45')
    ax[2,0].set_ylabel('SSP5-85')


def map_plot(cube, ax, levels, cmap):
    lats=cube.coord('latitude').points
    lons=cube.coord('longitude').points
    c = ax.contourf(lons,lats,cube.data, levels = levels, cmap = cmap, transform=ccrs.PlateCarree(), extend='both')
    ax.add_feature(cf.BORDERS)
    ax.set_xlim([79, 90])
    ax.set_ylim([23, 33])
    #gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    #ax.colorbar(c)
    return c

def ts_plot(ts, lab_text, title_text):
    ''' Plot relative to first decade'''
    fig = plt.figure()
    for each, label_text in zip(ts, lab_text):
        icc.add_year(each, each.coord('time'))
        each.data = each.data - np.mean(each.data[:10])
        plt.plot(each.coord('year').points, each.data, label=label_text)
    plt.axvspan(2041, 2060, color='red', alpha=0.3)
    plt.axvspan(2081, 2100, color='red', alpha=0.3)
    plt.xlabel('Year')
    plt.ylabel('Temp, *C')
    plt.title(title_text)
    plt.legend()
    return

def plot_box_on_map(bdry):
    plt.plot([bdry[3]. bdry[2]], [bdry[1], bdry[1]],'k')
    plt.plot([bdry[3]. bdry[2]], [bdry[0], bdry[0]],'k')
    plt.plot([bdry[3]. bdry[3]], [bdry[1], bdry[0]],'k')
    plt.plot([bdry[2]. bdry[2]], [bdry[1], bdry[0]],'k')
    return

def clear_tmp():
    files = glob.glob('tmp/*.nc')
    for f in files:
        os.remove(f)
    return

def timeslices(infile, seas, label, region):
    cube = seasonal_timeseries(infile, seas, label)
    icc.add_year(cube, cube.coord('time'))
    cube1 = cube.extract(iris.Constraint(year=lambda cell:2021 <= cell <= 2040))
    cube2 = cube.extract(iris.Constraint(year=lambda cell:2041 <= cell <= 2060))
    cube3 = cube.extract(iris.Constraint(year=lambda cell:2081 <= cell <= 2100))
    cube1 = extract_box(cube1[0,:,:], region)-273.15
    cube2 = extract_box(cube2[0,:,:], region)-273.15
    cube3 = extract_box(cube3[0,:,:], region)-273.15
    return cube1, cube2, cube3

def timeslices_rel(infile, seas, label, region):
    cube = seasonal_timeseries(infile, seas, label)
    icc.add_year(cube, cube.coord('time'))
    cube1 = cube.extract(iris.Constraint(year=lambda cell:2021 <= cell <= 2030))
    cube1 = extract_box(cube1[0,:,:], region)-273.15
    return cube1

def timeseries_plot(region, seas, title):
    clear_tmp()
    cube = seasonal_timeseries(infile585, seas, '585')
    new_cube = extract_box(cube, region)
    ts585 = new_cube.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)-273.15
    cube = seasonal_timeseries(infile245, seas, '245')
    new_cube = extract_box(cube, region)
    ts245 = new_cube.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)-273.15
    cube = seasonal_timeseries(infile126, seas, '126')
    new_cube = extract_box(cube, region)
    ts126 = new_cube.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)-273.15
    ts = [ts585, ts245, ts126]
    labs = ['ssp5-85', 'ssp2-45', 'ssp1-26']
    ts_plot(ts, labs, title+seas)
    return ts

def plots(reg_main, seas, reg_list):
    ## Timeslice map
    clear_tmp()
    cube1, cube2, cube3 = timeslices(infile126, seas, '126', reg_main)
    cube4, cube5, cube6 = timeslices(infile245, seas, '245', reg_main)
    cube7, cube8, cube9 = timeslices(infile585, seas, '585', reg_main)
    cubes = [cube1, cube2, cube3, cube4, cube5, cube6, cube7, cube8, cube9]
    multiple_maps(cubes, np.linspace(0, 45, 20), plt.cm.get_cmap('Reds')) # linspace bit defines color range
    # Rel to 2021-2040
    clear_tmp()
    cube1, cube2, cube3 = timeslices(infile126, seas, '126', reg_main)
    cube4, cube5, cube6 = timeslices(infile245, seas, '245', reg_main)
    cube7, cube8, cube9 = timeslices(infile585, seas, '585', reg_main)
    cubes = [cube1, cube2, cube3, cube4, cube5, cube6, cube7, cube8, cube9]
    rel_cube = timeslices_rel(infile126, seas, '126', reg_main)
    new_cubes = []
    for each in cubes:
        new_cubes.append(each-rel_cube)
    multiple_maps(new_cubes, np.linspace(-5, 5, 20), plt.cm.get_cmap('RdBu_r')) 
    ## Timeseries plot
    ts = timeseries_plot(reg_main, seas, 'Nepal, ')
    ts_west = timeseries_plot(reg_list[0], seas, 'Nepal west, ')
    ts_east = timeseries_plot(reg_list[1], seas, 'Nepal east, ')



##########################
## IPCC
# Near term: 2021-2040
# Medium term: 2041-2060
# Long term: 2081-2100
# Baseline 1981-2010

# HadGEM as an example
infile585 = '/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/MOHC/HadGEM3-GC31-LL/ssp585/r1i1p1f3/Amon/tas/gn/latest/*'
infile126 = '/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/MOHC/HadGEM3-GC31-LL/ssp126/r1i1p1f3/Amon/tas/gn/latest/*'
infile245 = '/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/MOHC/HadGEM3-GC31-LL/ssp245/r1i1p1f3/Amon/tas/gn/latest/*'

reg1 = [23, 33, 78, 91] # Nepal
reg_west = [23, 33, 78, 85] # Nepal west
reg_east = [23, 33, 85, 91] # Nepal east
plots(reg1, 'DJF', [reg_west, reg_east])
plots(reg1, 'JJA', [reg_west, reg_east])



# have got day/tasmax   and day/pr






### Seasonal precipitation


def seasonal_timeseries(infile, seas, label):
    ''' seas = 'JJA'  years = 'XXXX' '''
    ''' '''
    # Merge to one file 2015-2100
    cdo_cmd = ['cdo', 'mergetime', infile, 'tmp/ssp'+label+'_p1.nc']
    ret = subprocess.call(cdo_cmd)
    # Extract JJA months
    cdo_seas_p1 = ['cdo', 'selseas,'+seas, 'tmp/ssp'+label+'_p1.nc', 'tmp/ssp'+label+'_p2.nc']
    ret = subprocess.call(cdo_seas_p1)
    # JJA annual values
    cdo_seas_p2 = ['cdo', 'seassum', 'tmp/ssp'+label+'_p2.nc', 'tmp/ssp_'+label+'.nc']
    ret = subprocess.call(cdo_seas_p2)
    return iris.load('tmp/ssp_'+label+'.nc')[0]

def multiple_maps(cubes, levels, cmap):
    fig, ax = plt.subplots(nrows=3, ncols=3, subplot_kw={'projection': ccrs.PlateCarree()})
    ax_list = [ax[0,0], ax[0,1], ax[0,2], ax[1,0], ax[1,1], ax[1,2], ax[2,0], ax[2,1], ax[2,2]]
    for i, each in enumerate(cubes):
        c = map_plot(each, ax_list[i], levels, cmap)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(c, cax=cbar_ax)
    ax[0,0].set_title('2021-2040')
    ax[0,1].set_title('2041-2060')
    ax[0,2].set_title('2081-2100')
    ax[0,0].set_ylabel('SSP1-26')
    ax[1,0].set_ylabel('SSP2-45')
    ax[2,0].set_ylabel('SSP5-85')


def map_plot(cube, ax, levels, cmap):
    lats=cube.coord('latitude').points
    lons=cube.coord('longitude').points
    c = ax.contourf(lons,lats,cube.data, levels = levels, cmap = cmap, transform=ccrs.PlateCarree(), extend='both')
    ax.add_feature(cf.BORDERS)
    ax.set_xlim([79, 90])
    ax.set_ylim([23, 33])
    #gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    #ax.colorbar(c)
    return c

def ts_plot(ts, lab_text, title_text):
    ''' Plot relative to first decade'''
    fig = plt.figure()
    for each, label_text in zip(ts, lab_text):
        icc.add_year(each, each.coord('time'))
        each.data = each.data - np.mean(each.data[:10])
        plt.plot(each.coord('year').points, each.data, label=label_text)
    plt.axvspan(2041, 2060, color='red', alpha=0.3)
    plt.axvspan(2081, 2100, color='red', alpha=0.3)
    plt.xlabel('Year')
    plt.ylabel('Precipitation')
    plt.title(title_text)
    plt.legend()
    return

def timeslices(infile, seas, label, region):
    cube = seasonal_timeseries(infile, seas, label)
    icc.add_year(cube, cube.coord('time'))
    cube1 = cube.extract(iris.Constraint(year=lambda cell:2021 <= cell <= 2040))
    cube2 = cube.extract(iris.Constraint(year=lambda cell:2041 <= cell <= 2060))
    cube3 = cube.extract(iris.Constraint(year=lambda cell:2081 <= cell <= 2100))
    cube1 = extract_box(cube1[0,:,:], region)
    cube2 = extract_box(cube2[0,:,:], region)
    cube3 = extract_box(cube3[0,:,:], region)
    return cube1, cube2, cube3

def timeslices_rel(infile, seas, label, region):
    cube = seasonal_timeseries(infile, seas, label)
    icc.add_year(cube, cube.coord('time'))
    cube1 = cube.extract(iris.Constraint(year=lambda cell:2021 <= cell <= 2030))
    cube1 = extract_box(cube1[0,:,:], region)
    return cube1

def timeseries_plot(region, seas, title):
    clear_tmp()
    cube = seasonal_timeseries(infile585, seas, '585')
    new_cube = extract_box(cube, region)
    ts585 = new_cube.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)
    cube = seasonal_timeseries(infile245, seas, '245')
    new_cube = extract_box(cube, region)
    ts245 = new_cube.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)
    cube = seasonal_timeseries(infile126, seas, '126')
    new_cube = extract_box(cube, region)
    ts126 = new_cube.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)
    ts = [ts585, ts245, ts126]
    labs = ['ssp5-85', 'ssp2-45', 'ssp1-26']
    ts_plot(ts, labs, title+seas)
    return ts

def plots(reg_main, seas, reg_list):
    ## Timeslice map
    clear_tmp()
    cube1, cube2, cube3 = timeslices(infile126, seas, '126', reg_main)
    cube4, cube5, cube6 = timeslices(infile245, seas, '245', reg_main)
    cube7, cube8, cube9 = timeslices(infile585, seas, '585', reg_main)
    cubes = [cube1, cube2, cube3, cube4, cube5, cube6, cube7, cube8, cube9]
    multiple_maps(cubes, np.linspace(0, 45, 20), plt.cm.get_cmap('Reds')) # linspace bit defines color range
    # Rel to 2021-2040
    clear_tmp()
    cube1, cube2, cube3 = timeslices(infile126, seas, '126', reg_main)
    cube4, cube5, cube6 = timeslices(infile245, seas, '245', reg_main)
    cube7, cube8, cube9 = timeslices(infile585, seas, '585', reg_main)
    cubes = [cube1, cube2, cube3, cube4, cube5, cube6, cube7, cube8, cube9]
    rel_cube = timeslices_rel(infile126, seas, '126', reg_main)
    new_cubes = []
    for each in cubes:
        new_cubes.append(each-rel_cube)
    multiple_maps(new_cubes, np.linspace(-5, 5, 20), plt.cm.get_cmap('RdBu_r')) 
    ## Timeseries plot
    ts = timeseries_plot(reg_main, seas, 'Nepal, ')
    ts_west = timeseries_plot(reg_list[0], seas, 'Nepal west, ')
    ts_east = timeseries_plot(reg_list[1], seas, 'Nepal east, ')



infile585 = '/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/MOHC/HadGEM3-GC31-LL/ssp585/r1i1p1f3/Amon/ps/gn/latest/*'
infile126 = '/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/MOHC/HadGEM3-GC31-LL/ssp126/r1i1p1f3/Amon/ps/gn/latest/*'
infile245 = '/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/MOHC/HadGEM3-GC31-LL/ssp245/r1i1p1f3/Amon/ps/gn/latest/*'

reg1= [23, 33, 78, 91] # Nepal
ts = timeseries_plot(reg1, 'JJA', 'Nepal, ')
plt.ylim([7000, 10000]) # DJF
plt.ylim([-300, 1500]) # JJA

reg_west = [23, 33, 78, 85] # Nepal west
reg_east = [23, 33, 85, 91] # Nepal east
plots(reg1, 'DJF', [reg_west, reg_east])
plots(reg1, 'JJA', [reg_west, reg_east])
