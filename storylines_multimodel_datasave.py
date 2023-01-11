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


# All avail on BP
model_centre = ['AS-RCEC/TaiESM1', 'CCCR-IITM/IITM-ESM', 'E3SM-Project/E3SM-1-1', 'MIROC/MIROC6', 'NCC/NorESM2-LM', 'AWI/AWI-CM-1-1-MR',  'CNRM-CERFACS/CNRM-ESM2-1', 'EC-Earth-Consortium/EC-Earth3','MOHC/HadGEM3-GC31-LL', 'NIMS-KMA/KACE-1-0-G', 'BCC/BCC-CSM2-MR', 'CSIRO/ACCESS-ESM1-5', 'FIO-QLNM/FIO-ESM-2-0', 'MPI-M/MPI-ESM1-2-LR', 'NOAA-GFDL/GFDL-ESM4', 'CAMS/CAMS-CSM1-0', 'MRI/MRI-ESM2-0','NUIST/NESM3','CAS/FGOALS-f3-L', 'INM/INM-CM5-0', 'NASA-GISS/GISS-E2-1-G','THU/CIESM','CCCma/CanESM5', 'DWD/MPI-ESM1-2-HR','IPSL/IPSL-CM6A-LR', 'NCAR/CESM2']

var_measure = {'tas': 'seasmean', 'pr':'seassum', 'tasmax':'seasmax', 'prmax':'seasmax'}
var_path = {'tas': '/Amon/tas', 'pr':'/Amon/pr', 'tasmax':'/day/tasmax', 'prmax':'/day/pr'}
reg = {'Nepal':[23, 33, 78, 91], 'Nepal_west':[23, 33, 78, 85],'Nepal_east':[23, 33, 85, 91], 'Ethiopia':[3, 17, 30, 50], 'Ethiopia_lo':[3, 17, 40, 50], 'Ethiopia_hi':[3, 17, 30, 40], 'Uganda':[-4, 6, 28, 35], 'Uganda_north':[1, 6, 28, 35], 'Uganda_south':[-4, 1, 28, 35]} # S, N, W, E

### VARIABLES
region_list = ['Nepal', 'Nepal_east', 'Nepal_west', 'Ethiopia', 'Ethiopia_lo', 'Ethiopia_hi', 'Uganda', 'Uganda_north', 'Uganda_south']
scenario_list = ['ssp126', 'ssp245', 'ssp585']   
climate_var = 'tasmax' # 'tas', 'pr', 'tasmax', 'prmax'
seas = 'DJF'

''' 
N.B. need to change in script for seas
pr / prmax - mm/day
tas / tasmax - *C 
Currently both relative to first decade - biases in models prevent others. Could plot ensemble mean absolute.
'''


# want to add maps of ensemble mean too
# should be able to add this in after the cube stage - no need to rerun that bit
# first will need mean of cubes, may nto all be on same grid?
# then code for map already there further down - just two timescales and one scenario

def get_cube_list(scenario, climate_var, seas):
    clear_tmp()
    cube_list = []
    lab_val = 50000
    for centre in model_centre:
        print(centre)
        infile_path = '/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/'+centre+'/'+scenario+'/'
        if os.path.isdir(infile_path) == True:
            ens_list = os.listdir(infile_path)
            if len(ens_list) > 3: 
                ens_list = ens_list[:2]
            for each in ens_list:
                infile = infile_path+each+var_path[climate_var]+'/g*/latest/*'
                if glob.glob(infile) != []:
                    print(infile)
                    lab_val = lab_val+1
                    cube_list.append(seasonal_timeseries(infile, seas, str(lab_val), var_measure[climate_var]))
                else:
                    pass
    return cube_list


def get_region_ts(region, cube_list):
    # for each in cube list extract region
    ts_list = []
    for each in cube_list:
        new_cube = extract_box(each, reg[region])
        if climate_var in {'tas','tasmax'}:
            ts = new_cube.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)-273.15
        elif climate_var == 'pr' or 'prmax':
            ts = new_cube.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)*24*60*60*90
        elif climate_var == 'prmax':
            ts = new_cube.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)*24*60*60
        ts_list.append(ts)
    return ts_list

def plot_spagetti(ts_list, region, scenario, climate_var, seas):
    # plot relative to first decade, and as 10-year smoothed'''
    fig = plt.figure()
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
                each.data = each.data - np.mean(each.data[:10])
                # 10-year rolling
                dec = []
                for i in np.arange((len(each.data)-10)):
                    dec.append(np.mean(each.data[i:i+10]))
                # +5 to make centre of decade plotted
                plt.plot(each.coord('year').points[:-10]+5, dec, color='k', alpha=0.2)
                dec_data.append(dec)
                dec_year.append(each.coord('year').points[:-10]+5)
            except AttributeError:
                pass
    # plot ensemble mean
    new_data = []
    for n, each in enumerate(dec_data):
        for i in range(dec_year[n][0]-2020):
            each.insert(0, np.nan)
        for i in range(2096-dec_year[n][-1]): # append nan this many times
            each = np.append(each, np.nan)
        new_data.append(each)
    plt.plot(np.arange(2020, 2097, 1), np.nanmean(new_data, axis = 0), color='r')    
    # format plot
    plt.xlim([2020, 2095])
    plt.axvspan(2045, 2055, color='red', alpha=0.3)
    plt.axvspan(2085, 2095, color='red', alpha=0.3)
    plt.xlabel('Year')
    plt.ylabel('Temp, *C')
    plt.title(region+'_'+scenario+'_'+climate_var+'_'+seas)
    plt.grid()
    plt.savefig(region+'_'+scenario+'_'+climate_var+'_'+seas+'.png')
    return 






def seasonal_timeseries(infile, seas,label, measure):
    ''' seas = 'JJA'  years = 'XXXX' '''
    ''' '''
    # Merge to one file 2015-2100
    try:
        cdo_cmd = ['cdo', 'mergetime', infile, 'tmp/ssp'+label+'_p1.nc']
        ret = subprocess.call(cdo_cmd)
        cdo_seas_p1 = ['cdo', 'selseas,'+seas, 'tmp/ssp'+label+'_p1.nc', 'tmp/ssp'+label+'_p2.nc']
        ret = subprocess.call(cdo_seas_p1)
        cdo_seas_p2 = ['cdo', measure, 'tmp/ssp'+label+'_p2.nc', 'tmp/ssp_'+label+'.nc']
        ret = subprocess.call(cdo_seas_p2) 
    except OSError:
        pass
    x = iris.load('tmp/ssp_'+label+'.nc')[0]
    os.remove('tmp/ssp'+label+'_p1.nc')
    os.remove('tmp/ssp'+label+'_p2.nc')
    return x
    

def extract_box(cube_in, bdry):
    ' bdry = list, [north, south, east, west] '
    const_lat = iris.Constraint(latitude = lambda cell:bdry[0] < cell < bdry[1])
    const_lon = iris.Constraint(longitude = lambda cell:bdry[2] < cell < bdry[3])
    new_cube = cube_in.extract(const_lat & const_lon)
    return new_cube

def clear_tmp():
    files = glob.glob('tmp/*.nc')
    for f in files:
        os.remove(f)
    return




seas = 'JJA'
climate_var = 'tasmax'
scenario = 'ssp245'
x = get_cube_list(scenario, climate_var, seas)
for region in region_list:
    nep_list = get_region_ts(region, x)
    plot_spagetti(nep_list, region, scenario, climate_var, seas)


climate_var = 'prmax'
scenario = 'ssp245'
for seas in ['JJA', 'DJF', 'SON']:
    x = get_cube_list(scenario, climate_var, seas)
    for region in region_list[:3]:
        nep_list = get_region_ts(region, x)
        plot_spagetti(nep_list, region, scenario, climate_var, seas)

for seas in ['MAM', 'SON']:
    x = get_cube_list(scenario, climate_var, seas)
    for region in region_list[6:]:
        nep_list = get_region_ts(region, x)
        plot_spagetti(nep_list, region, scenario, climate_var, seas)







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

