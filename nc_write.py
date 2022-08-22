# -*- coding: utf-8 -*-
import os
from pathlib import Path
import time as time_mod
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from conversions import VPD2RH
from conversions import convert_pr
from conversions import convert_ps
from conversions import convert_vpd
from conversions import convert_ta

FLUXNET_REF = """Pastorello, G., Trotta, C., Canfora, E. et al. The FLUXNET2015 dataset and the ONEFlux processing pipeline for eddy covariance data. Sci Data 7, 225 (2020). https://doi.org/10.1038/s41597-020-0534-3"""

# Sites metadata
SITES_COORDINATES = {'SITE': ['Latitude', 'Longitude', 'name', 'filepath'],
                     'Dav': [np.array([46.81], 'f4'), np.array([9.85], 'f4'), 'FLX-Dav',
                             './FLX_CH-Dav_FLUXNET2015_FULLSET_1997-2014_1-4/FLX_CH-Dav_FLUXNET2015_FULLSET_DD_1997-2014_1-4.csv'],
                     'Tha': [np.array([50.96], 'f4'), np.array([13.56], 'f4'), 'FLX-Tha',
                             './FLX_DE-Tha_FLUXNET2015_FULLSET_1996-2014_1-4/FLX_DE-Tha_FLUXNET2015_FULLSET_DD_1996-2014_1-4.csv']}

# Variables METAdata (Fprcing)
FLUXNET_FULLSET_VARS = {'tas':  ["TA_F_MDS_DAY", "Air temperature, gapfilled using MDS method", "K", 'air_temperature'], # FLUXNET celsius
                        'rsds': ["SW_IN_F_MDS", "Shortwave radiation, incoming, gapfilled using MDS", 'W m-2', "surface_downwelling_shortwave_flux_in_air"],
                        'vpd':  ['VPD_F_MDS', "Vapor Pressure Deficit gapfilled using MDS method", "kPa", "vpd"], # FLUXNET hPa
                        'ps': ["PA_F","Atmospheric Pressure","Pa", "surface_air_pressure"], # FLUXNET kPa
                        'pr': ["P_F","Precipitation","kg m-2 s-1", "precipitation_flux"], # FLUXNET mm/day
                        'wind': ["WS_F","Wind Speed","m s-1", "wind_speed"], # Aparently the same 
                        'PPFD': ["PPFD_IN","Photosynthetic photon flux density, incoming","µmol m-2 s-1", "PPFD"],
                        'co2': ["CO2_F_MDS","CO2 mole fraction, gapfilled with MDS", "µmol mol-1", "CONVERT TO TEXT\n<year> <value>\n"],
                        'hurs': ["RH_F","Relative humidity, range 0-100", "%", 'relative_humidity']}

# TODO Collect imformation about the reference variables (e.g. NEE, LH) for the benchmark

def save_gridlist(site="SITE"):
    lat = SITES_COORDINATES[site][0][0]
    lon = SITES_COORDINATES[site][1][0]
    fname = SITES_COORDINATES[site][2]
    
    with open(f"{fname}.txt", 'w', encoding="utf-8") as fh:
        fh.write(f"{lon}\t{lat}\t{fname}\r\n")

def create_dataset(fname = None, 
                   arr=None, 
                   var=None,
                   unit=None, 
                   descr=None, 
                   time=None, 
                   la=None, 
                   lo=None, 
                   site_data=True, 
                   set="FULLSET",
                   reference=None):
    
    """Create a CF (not so) compliant nc4 file for site(s) data TODO document """
    
    if fname is not None:  
        dset = Dataset(os.path.join(Path('./'), fname + ".nc"),mode="w", format="NETCDF4")
    else:
        assert var is not None, "Need to set fname or var"
        dset = Dataset(os.path.join(Path('./'), FLUXNET_FULLSET_VARS[var][0] + '.nc'), mode="w", format="NETCDF4")
        
    if site_data:
        # Create temporal dimension
        time_dim = time['data']
        time_unit = time['units']
        calendar = time['calendar']
        
        lats = la
        lons = lo
        
        # Create netCDF dimensions
        dset.createDimension("time",None)
        dset.createDimension("lon",size=lons.size)
        dset.createDimension("lat",size=lats.size)
        
  
        # Data description
        dset.description = f"FLUXNET 2015 DATA - {set} {FLUXNET_FULLSET_VARS[var][0]}"
        dset.source = f'Forcing data for DVM - {descr}'
        dset.history= f'Created: {time_mod.ctime(time_mod.time())}'
        if reference is not None:
            dset.reference = reference

        # Create netCDF variables

        X  = dset.createVariable("lon", 'f4' , ("lon",))
        Y =  dset.createVariable("lat", 'f4' , ("lat",))
        T  = dset.createVariable("time", 'i4', ("time",))
        D  = dset.createVariable(var, 'f8', ("time", "lon", "lat"))

        T[:] = time_dim.__array__()
        T.units    = time_unit
        T.calendar = calendar
        T.standard_name = "time"
        T.axis = 'T'
        
        X[:] = lons
        X.units    = "degrees_east"
        X.long_name = 'longitude'
        X.standard_name = 'longitude'
        X.axis = 'X'

        Y[:] = lats
        Y.units    = "degrees_north"
        Y.long_name = 'latitude'
        Y.standard_name = 'latitude'
        Y.axis = 'Y' 
               
        D[:, 0, 0] = np.copy(arr[:], order="C")
        D.units = unit
        # D.fluxnet_name = FLUXNET_FULLSET_VARS[var][0]
        D.standard_name = FLUXNET_FULLSET_VARS[var][3]
        D.coordinates = u"time lon lat"
        # D.coordinates = 'lon lat'
        
        dset.close()
        
    else:
        # TODO implement layered data
        pass

def single_timeseries(fname = None, 
                      arr=None, 
                      var=None,
                      site=None,
                      unit=None, 
                      descr=None, 
                      time=None, 
                      la=None, 
                      lo=None, 
                      site_data=True, 
                      set="FULLSET",
                      reference=None):
    
    """Create a CF compliant nc4 file for site(s) data TODO document
    Single Timeseries CF conventions (page 163)"""
    
    
    if fname is not None:  
        dset = Dataset(os.path.join(Path('./'), fname + ".nc"),mode="w", format="NETCDF4")
    else:
        assert var is not None, "Need to set fname or var"
        dset = Dataset(os.path.join(Path('./'), FLUXNET_FULLSET_VARS[var][0] + '.nc'), mode="w", format="NETCDF4")
        
    if site_data:
        # Create temporal dimension
        time_dim = time['data']
        time_unit = time['units']
        calendar = time['calendar']
        
        lats = la
        lons = lo
        
        # Create netCDF dimensions
        dset.createDimension("time",size=arr.size)
        # dset.createDimension("station_name",size=7)
        #dset.createDimension("lon",size=lons.size)
        #dset.createDimension("lat",size=lats.size)
        
  
        # Data description
        dset.description = f"FLUXNET 2015 DATA - {set} {FLUXNET_FULLSET_VARS[var][0]}"
        dset.source = f'Forcing data for DVM - {descr}'
        dset.history= f'Created: {time_mod.ctime(time_mod.time())}'
        if reference is not None:
            dset.reference = reference
        dset.featureType = "timeSeries"

        # Create netCDF variables
        X  = dset.createVariable("lon", 'f4')
        Y =  dset.createVariable("lat", 'f4')
        SN = dset.createVariable("station_name", '<U7')
        T  = dset.createVariable("time", 'i4', ("time",))
        
        D  = dset.createVariable(var, 'f8', ("time",))

        T[...] = time_dim.__array__()
        T.units    = time_unit
        T.calendar = calendar
        T.standard_name = "time"
        T.axis = 'T'
        
        X[...] = lons[0]
        X.units    = "degrees_east"
        X.long_name = 'station_longitude'
        X.standard_name = 'longitude'

        Y[...] = lats[0]
        Y.units    = "degrees_north"
        Y.long_name = 'station_latitude'
        Y.standard_name = 'latitude'
        
        SN[...] = SITES_COORDINATES[site][2]
        SN.long_name = "station name"
        SN.cf_role = "timeseries_id"
        
                      
        D[...] = np.copy(arr[:], order="C")
        D.units = unit
        # D.fluxnet_name = FLUXNET_FULLSET_VARS[var][0]
        D.standard_name = FLUXNET_FULLSET_VARS[var][3]
        D.coordinates = u"time lat lon station_name"
        # D.coordinates = 'lon lat'
        
        dset.close()
        
    else:
        # TODO implement layered data
        pass

def write_site_nc(VAR, SITE):
    """write FLUXNET2015 FULLSET Variable to a netCDF4 file """
    fluxnet_data = pd.read_csv(SITES_COORDINATES[SITE][3])
    last = fluxnet_data['TIMESTAMP'][:].size - 1
    start = str(fluxnet_data['TIMESTAMP'][0])
    end = str(fluxnet_data['TIMESTAMP'][last])

    idx = pd.date_range(start, end, freq='D')

    # TODO implement standard calendar "1850-01-01T00:00:00" here
    time_data = np.arange(idx.size, dtype='i4')

    # mid_day = pd.offsets.Hour() * 12
    #day_init, hour_init = (idx[0] + mid_day).isoformat().split("T")
    day_init, hour_init = idx[0].isoformat().split("T")
    time_units = "days since %s %s" % (str(day_init), str(hour_init))
    calendar = 'proleptic_gregorian'

    descr = FLUXNET_FULLSET_VARS[VAR][1]
    units = FLUXNET_FULLSET_VARS[VAR][2]

    time_dict = {'data': time_data,
                'units': time_units,
                'calendar': calendar}

    success = False
    try:
        if VAR == 'pr':
            arr = convert_pr(fluxnet_data[FLUXNET_FULLSET_VARS[VAR][0]].__array__())
        elif VAR == 'ps':
            arr = convert_ps(fluxnet_data[FLUXNET_FULLSET_VARS[VAR][0]].__array__())
        elif VAR == 'vpd':
            arr = convert_vpd(fluxnet_data[FLUXNET_FULLSET_VARS[VAR][0]].__array__())
        elif VAR == 'tas':
            arr = convert_ta(fluxnet_data[FLUXNET_FULLSET_VARS[VAR][0]].__array__())
        else:
            arr = fluxnet_data[FLUXNET_FULLSET_VARS[VAR][0]].__array__()
        success = True
    except:
        if VAR == 'hurs':
            try:
                vpd = fluxnet_data[FLUXNET_FULLSET_VARS['vpd'][0]].__array__()
                tair = fluxnet_data[FLUXNET_FULLSET_VARS['tas'][0]].__array__()
                arr = 100.0 * VPD2RH(tair, vpd * 0.1)
                success = True
            except:
                success = False
        pass
    
    if success:
        assert arr.size == time_data.size, "arr dif size"
        la = SITES_COORDINATES[SITE][0]
        lo = SITES_COORDINATES[SITE][1]
        
        fname00 = f"{VAR}_{SITES_COORDINATES[SITE][2]}"
        f3 = 'out0'
        
        # single_timeseries(fname=fname00, arr=arr, var=VAR, site=SITE,unit=units, 
        #                descr=descr, time=time_dict, la=la, lo=lo,
        #                reference=FLUXNET_REF)
        

        create_dataset(fname=fname00, arr=arr, var=VAR, unit=units, 
                       descr=descr, time=time_dict, la=la, lo=lo,
                       reference=FLUXNET_REF)

        # Change calendar and etc.
        # os.system(f"cdo settaxis,{day_init},{hour_init},days TEST_FLUXNET_NETCDF.nc4 out.nc4")
        # os.system("cdo delete,month=5,7,8,10,12,day=31 out0.nc4 out1.nc4")
        
        # os.system("cdo delete,month=2,day=29 out0.nc4 out5.nc4")
        # os.system(f"cdo setcalendar,365_day out5.nc4  out6.nc4")
        # os.system(f"cdo settaxis,{day_init},{hour_init},day out6.nc4 out7.nc4")
        # os.system(f"cdo settbounds,day out7.nc4 {VAR}_{SITES_COORDINATES[SITE][2]}.nc4")
        # os.system("rm -rf out*")
        
        # Write co2 as annual text files
        # if VAR == 'co2':
        #     print(f"Processing CO2 {SITE}")
        #     os.system(f"cdo yearmean {VAR}_{SITES_COORDINATES[SITE][2]}.nc4 {VAR}_{SITES_COORDINATES[SITE][2]}_yearAvg.nc4")
        #     with Dataset(f"{VAR}_{SITES_COORDINATES[SITE][2]}_yearAvg.nc4", 'r') as fh:
        #         tindex = fh.variables['time'][:]
        #         calendar_aux = fh.variables['time'].calendar
        #         tunits_aux = fh.variables['time'].units
        #         co2Data = fh.variables['co2'][:]
        #         year0 = cftime.num2date(tindex[0], tunits_aux, calendar_aux)
        #         yearF = cftime.num2date(tindex[-1], tunits_aux, calendar_aux)
            
        #     lines = []
        #     for i, x in enumerate(range(year0.year, yearF.year + 1)):
        #         line_to_append = f"{x} {round(co2Data[i][0][0], 2)}\n" 
        #         lines.append(line_to_append)
        #         print(f"Created CO2 for the year {x} = {round(co2Data[i][0][0], 2)}")
            
            
        #     with open(f"co2_FLX_{SITE}.txt", 'w') as fh:
        #         fh.writelines(lines)
    else:
        print(f"VAR NOT FOUND: {FLUXNET_FULLSET_VARS[VAR][0]}")