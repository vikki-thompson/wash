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
reg = {'Nepal': [23, 33, 78, 91], 'Nepal_west':[23, 33, 78, 85],'Nepal_east':[23, 33, 85, 91], 'Ethiopia':[3, 17, 30, 50], 'Ethiopia_lo':[3, 17, 40, 50], 'Ethiopia_hi':[3, 17, 30, 40], 'Uganda':[-4, 6, 28, 35], 'Uganda_north':[0, 6, 28, 35], 'Uganda_south':[-4, 1, 28, 35], 'Senegal':[11, 18, 341, 349]} # S, N, W, E


def get_cube_list(scenario, climate_var, seas):
    clear_tmp()
    cube_list = []
    lab_val = 50000
    for centre in model_centre:
        print(centre)
        #if centre in {'AWI/AWI-CM-1-1-MR','EC-Earth-Consortium/EC-Earth3','NIMS-KMA/KACE-1-0-G', 'CSIRO/ACCESS-ESM1-5','MRI/MRI-ESM2-0', 'IPSL/IPSL-CM6A-LR', 'MPI-M/MPI-ESM1-2-LR'} and climate_var in {'prmax', 'tasmax'}: 
        #    pass
        #else:
        infile_path = '/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/'+centre+'/'+scenario+'/'
        if os.path.isdir(infile_path) == True:
            ens_list = os.listdir(infile_path)
            if len(ens_list) > 3: 
                ens_list = ens_list[:2]
            for each in ens_list:
                #infile = infile_path+each+var_path[climate_var]+'/g*/latest/*'
                infile = infile_path+each+'/day/pr'+'/g*/latest/*'
                lab_val = lab_val+1
                if glob.glob(infile) != []:
                    print(infile)
                    cube_list.append(seasonal_timeseries(infile, seas, str(lab_val), var_measure[climate_var]))
                else:
                    pass
    return cube_list


def get_cube_list_2day(scenario, climate_var, seas):
    clear_tmp()
    cube_list = []
    lab_val = 50032
    for centre in model_centre[13:]:  
        print(centre)
        infile_path = '/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/'+centre+'/'+scenario+'/'
        if os.path.isdir(infile_path) == True:
            ens_list = os.listdir(infile_path)
            print(ens_list)
            if len(ens_list) > 3: 
                ens_list = ens_list[:2]
            for each in ens_list:
                #infile = infile_path+each+var_path[climate_var]+'/g*/latest/*'
                infile = infile_path+each+'/day/pr/g*/latest/*'
                lab_val = lab_val+1
                print(lab_val)
                if glob.glob(infile) != []:
                    print(infile)
                    cube_list.append(seasonal_timeseries_2day(infile, seas, str(lab_val)))
                    iris.save(cube_list, '/user/work/hh21501/wash/pr2_'+str(lab_val)+'.nc')
                else:
                    pass
    return cube_list


'''
pr2_list = glob.glob('/user/work/hh21501/wash/pr2_5*.nc')
new_list = []
for each in pr2_list[:1]:
    x=iris.load(each)
    for c in x:
        if np.ndim(c) > 1:
            new_list.append(c)
        else: pass

iris.save(x, '/user/work/hh21501/wash/pr2_ssp242_jja.nc')
'''



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



def seasonal_timeseries(infile, seas, label, measure):
    ''' seas = 'JJA'  years = 'XXXX' '''
    ''' '''
    # Merge to one file 2015-2100
    try:
        cdo_cmd = ['cdo', 'mergetime', infile, '/user/work/hh21501/wash/tmp/ssp'+label+'_p1.nc']
        ret = subprocess.call(cdo_cmd)
        cdo_seas_p1 = ['cdo', 'selseas,'+seas, '/user/work/hh21501/wash/tmp/ssp'+label+'_p1.nc', '/user/work/hh21501/wash/tmp/ssp'+label+'_p2.nc']
        ret = subprocess.call(cdo_seas_p1)
        cdo_seas_p2 = ['cdo', measure, '/user/work/hh21501/wash/tmp/ssp'+label+'_p2.nc', '/user/work/hh21501/wash/tmp/ssp_'+label+'.nc']
        ret = subprocess.call(cdo_seas_p2) 
    except OSError:
        pass
    x = iris.load('/user/work/hh21501/wash/tmp/ssp_'+label+'.nc')[0].copy()
    #clear_tmp()
    return x


def seasonal_timeseries_5day(infile, seas, label):
    ''' seas = 'JJA'  years = 'XXXX' '''
    ''' '''
    # Merge to one file 2015-2100
    try:
        cdo_cmd = ['cdo', 'mergetime', infile, '/user/work/hh21501/wash/tmp/ssp'+label+'_p1.nc']
        ret = subprocess.call(cdo_cmd)
        cdo_seas_p1 = ['cdo', 'selseas,'+seas, '/user/work/hh21501/wash/tmp/ssp'+label+'_p1.nc', '/user/work/hh21501/wash/tmp/ssp'+label+'_p2.nc']
        ret = subprocess.call(cdo_seas_p1)
        x = iris.load('/user/work/hh21501/wash/tmp/ssp'+label+'_p2.nc')[0].copy()
        iris.coord_categorisation.add_year(x, 'time')
        sample = np.arange(x.coord('year').points[0], x.coord('year').points[-1]+1)
        output = x.aggregated_by('year', iris.analysis.MEAN)
        for i, yr in enumerate(sample):
            print(i, yr)
            yr_x = x.extract(iris.Constraint(year=yr))
            for lat in np.arange(0, np.shape(yr_x)[1]):
                for lon in np.arange(0, np.shape(yr_x)[2]):
                    ts = yr_x.data[:,lat,lon]
                    ts5 = []
                    for each in np.arange(0,len(ts.data)-4):
                        ts5.append(np.sum(ts[each:each+5].data))
                    output.data[i,lat,lon] = np.max(ts5)
    except OSError:
        pass
    return output

    
def seasonal_timeseries_2day(infile, seas, label):
    ''' seas = 'JJA'  years = 'XXXX' '''
    ''' '''
    # Merge to one file 2015-2100
    try:
        cdo_cmd = ['cdo', 'mergetime', infile, '/user/work/hh21501/wash/tmp/ssp'+label+'_p1.nc']
        ret = subprocess.call(cdo_cmd)
        cdo_seas_p1 = ['cdo', 'selseas,'+seas, '/user/work/hh21501/wash/tmp/ssp'+label+'_p1.nc', '/user/work/hh21501/wash/tmp/ssp'+label+'_p2.nc']
        ret = subprocess.call(cdo_seas_p1)
        x = iris.load('/user/work/hh21501/wash/tmp/ssp'+label+'_p2.nc')[0].copy()
        iris.coord_categorisation.add_year(x, 'time')
        sample = np.arange(x.coord('year').points[0], x.coord('year').points[-1]+1)
        output = x.aggregated_by('year', iris.analysis.MEAN)
        for i, yr in enumerate(sample):
            print(i, yr)
            yr_x = x.extract(iris.Constraint(year=yr))
            for lat in np.arange(0, np.shape(yr_x)[1]):
                for lon in np.arange(0, np.shape(yr_x)[2]):
                    ts = yr_x.data[:,lat,lon]
                    ts5 = []
                    for each in np.arange(0,len(ts.data)-1):
                        ts5.append(np.sum(ts[each:each+2].data))
                    output.data[i,lat,lon] = np.max(ts5)
    except OSError:
        pass
    return output


def extract_box(cube_in, bdry):
    ' bdry = list, [north, south, east, west] '
    const_lat = iris.Constraint(latitude = lambda cell:bdry[0] < cell < bdry[1])
    const_lon = iris.Constraint(longitude = lambda cell:bdry[2] < cell < bdry[3])
    new_cube = cube_in.extract(const_lat & const_lon)
    return new_cube


def clear_tmp():
    files = glob.glob('/user/work/hh21501/wash/tmp/*.nc')
    for f in files:
        os.remove(f)
    return






### VARIABLES
region_list = ['Nepal', 'Nepal_west', 'Nepal_east']
scenario_list = ['ssp245']   
climate_var = 'tas' # 'tas', 'pr', 'tasmax', 'prmax'
seas = 'DJF'

''' 
N.B. need to chnage in script for seas
pr / prmax - mm/day
tas / tasmax - *C 
Currently both relative to first decade - biases in models prevent others. Could plot ensemble mean absolute.
'''

## Get data and save cubes (only run once)
#x = get_cube_list('ssp245', 'tas', 'DJF')
#iris.save(x, '/user/work/hh21501/wash/tas_ssp242_djf.nc')
#x = get_cube_list('ssp245', 'tas', 'JJA')
#iris.save(x, '/user/work/hh21501/wash/tas_ssp242_jja.nc')
x = get_cube_list('ssp245', 'pr', 'MAM')
iris.save(x, '/user/work/hh21501/wash/pr_ssp242_mam2.nc')
x = get_cube_list('ssp245', 'pr', 'OND')
iris.save(x, '/user/work/hh21501/wash/pr_ssp242_ond2.nc')

##x = get_cube_list('ssp245', 'tasmax', 'MAM')
##iris.save(x, '/user/work/hh21501/wash/tasmax_ssp242_mam.nc')

#x = get_cube_list('ssp245', 'prmax', 'JJA')
#iris.save(x, '/user/work/hh21501/wash/prmax_ssp242_jja.nc')

#x = get_cube_list_2day('ssp245', 'prmax', 'JJA')
#iris.save(x, '/user/work/hh21501/wash/pr2_ssp242_jja.nc')
