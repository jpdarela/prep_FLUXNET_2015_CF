# -*- coding: utf-8 -*-
import os
from pathlib import Path
import time as time_mod
from netCDF4 import Dataset
import cftime
import numpy as np
import pandas as pd
from conversions import convert_pr, convert_ps, convert_ta, convert_vpd, VPD2RH


FLUXNET_REF = """Pastorello, G., Trotta, C., Canfora, E. et al. The FLUXNET2015 dataset and the ONEFlux processing pipeline for eddy covariance data. Sci Data 7, 225 (2020). https://doi.org/10.1038/s41597-020-0534-3"""

# # Sites metadata
# ERA downscaled
SITES_COORDINATES = {'SITE': ['Latitude', 'Longitude', 'name', 'filepath'],
                     'Dav': [np.array([46.8153], 'f4'), np.array([9.8559], 'f4'), 'FLX-Dav',
                             './FLX_CH-Dav_FLUXNET2015_FULLSET_1997-2014_1-4/FLX_CH-Dav_FLUXNET2015_ERAI_DD_1989-2014_1-4.csv'],
                    
                     'Tha': [np.array([50.9626], 'f4'), np.array([13.5651], 'f4'), 'FLX-Tha',
                             './FLX_DE-Tha_FLUXNET2015_FULLSET_1996-2014_1-4/FLX_DE-Tha_FLUXNET2015_ERAI_DD_1989-2014_1-4.csv'],
                    
                     'Hai': [np.array([51.0792], 'f4'), np.array([10.4522], 'f4'), 'FLX-Hai',
                             './FLX_DE-Hai_FLUXNET2015_FULLSET_2000-2012_1-4/FLX_DE-Hai_FLUXNET2015_ERAI_DD_1989-2014_1-4.csv'],
                    
                     'Lnf': [np.array([51.3282], 'f4'), np.array([10.3678], 'f4'), 'FLX-Lnf',
                             './FLX_DE-Lnf_FLUXNET2015_FULLSET_2002-2012_1-4/FLX_DE-Lnf_FLUXNET2015_ERAI_DD_1989-2014_1-4.csv'],
                    
                     'Obe': [np.array([50.7867], 'f4'), np.array([13.7213], 'f4'), 'FLX-Obe',
                             './FLX_DE-Obe_FLUXNET2015_FULLSET_2008-2014_1-4/FLX_DE-Obe_FLUXNET2015_ERAI_DD_1989-2014_1-4.csv'],
                    
                     'Lae': [np.array([47.4783], 'f4'), np.array([8.3644], 'f4'), 'FLX-Lae',
                             './FLX_CH-Lae_FLUXNET2015_FULLSET_2004-2014_1-4/FLX_CH-Lae_FLUXNET2015_ERAI_DD_1989-2014_1-4.csv'],
                     
                     'BK1': [np.array([49.5021], 'f4'), np.array([18.5369], 'f4'), 'FLX-BK1',
                             './FLX_CZ-BK1_FLUXNET2015_FULLSET_2004-2014_2-4/FLX_CZ-BK1_FLUXNET2015_ERAI_DD_1989-2014_2-4.csv'],
                     
                     'Lkb': [np.array([49.0996], 'f4'), np.array([13.3047], 'f4'), 'FLX-Lkb',
                             './FLX_DE-Lkb_FLUXNET2015_FULLSET_2009-2013_1-4/FLX_DE-Lkb_FLUXNET2015_ERAI_DD_1989-2014_1-4.csv'],
                     
                     'Sor': [np.array([55.4859], 'f4'), np.array([11.6446], 'f4'), 'FLX-Sor',
                             './FLX_DK-Sor_FLUXNET2015_FULLSET_1996-2014_2-4/FLX_DK-Sor_FLUXNET2015_ERAI_DD_1989-2014_2-4.csv'],
                     
                     'Col': [np.array([41.8494], 'f4'), np.array([13.5881], 'f4'), 'FLX-Col',
                             './FLX_IT-Col_FLUXNET2015_FULLSET_1996-2014_1-4/FLX_IT-Col_FLUXNET2015_ERAI_DD_1989-2014_1-4.csv'],
                    
                     'Ren': [np.array([46.5869], 'f4'), np.array([11.4337], 'f4'), 'FLX-Ren',
                             './FLX_IT-Ren_FLUXNET2015_FULLSET_1998-2013_1-4/FLX_IT-Ren_FLUXNET2015_ERAI_DD_1989-2014_1-4.csv'],
                     
                     'Fyo': [np.array([56.4615], 'f4'), np.array([32.9221], 'f4'), 'FLX-Fyo',
                             './FLX_RU-Fyo_FLUXNET2015_FULLSET_1998-2014_2-4/FLX_RU-Fyo_FLUXNET2015_ERAI_DD_1989-2014_2-4.csv']}

# # Variables METAdata (Fprcing)
FLUXNET_FULLSET_VARS = {'tas':  ["TA_ERA", "Air temperature, gapfilled using MDS method", "K", 'air_temperature'], # FLUXNET celsius
                        'rsds': ["SW_IN_ERA", "Shortwave radiation, incoming, gapfilled using MDS", 'W m-2', "surface_downwelling_shortwave_flux_in_air"],
                        'vpd':  ['VPD_ERA', "Vapor Pressure Deficit gapfilled using MDS method", "kPa", "vpd"], # FLUXNET hPa
                        'ps': ["PA_ERA","Atmospheric Pressure","Pa", "surface_air_pressure"], # FLUXNET kPa
                        'pr': ["P_ERA","Precipitation","kg m-2 s-1", "precipitation_flux"], # FLUXNET mm/day
                        'wind': ["WS_ERA","Wind Speed","m s-1", "wind_speed"], # Aparently the same
                        'PPFD': ["PPFD_IN","Photosynthetic photon flux density, incoming","µmol m-2 s-1", "PPFD"],
                        'co2': ["CO2_F_MDS","CO2 mole fraction, gapfilled with MDS", "µmol mol-1", "CONVERT TO TEXT\n<year> <value>\n"],
                        'hurs': ["RH_F","Relative humidity, range 0-100", "%", 'relative_humidity']}

# Reference data
OBS_SITES = {'SITE': 'filepath',
             'Dav': './FLX_CH-Dav_FLUXNET2015_FULLSET_1997-2014_1-4/FLX_CH-Dav_FLUXNET2015_FULLSET_MM_1997-2014_1-4.csv',
             'Tha': './FLX_DE-Tha_FLUXNET2015_FULLSET_1996-2014_1-4/FLX_DE-Tha_FLUXNET2015_FULLSET_MM_1996-2014_1-4.csv',
             'Hai': './FLX_DE-Hai_FLUXNET2015_FULLSET_2000-2012_1-4/FLX_DE-Hai_FLUXNET2015_FULLSET_MM_2000-2012_1-4.csv',
             'Lnf': './FLX_DE-Lnf_FLUXNET2015_FULLSET_2002-2012_1-4/FLX_DE-Lnf_FLUXNET2015_FULLSET_MM_2002-2012_1-4.csv',
             'Obe': './FLX_DE-Obe_FLUXNET2015_FULLSET_2008-2014_1-4/FLX_DE-Obe_FLUXNET2015_FULLSET_MM_2008-2014_1-4.csv',
             'Lae': './FLX_CH-Lae_FLUXNET2015_FULLSET_2004-2014_1-4/FLX_CH-Lae_FLUXNET2015_FULLSET_MM_2004-2014_1-4.csv',
             'BK1': './FLX_CZ-BK1_FLUXNET2015_FULLSET_2004-2014_2-4/FLX_CZ-BK1_FLUXNET2015_FULLSET_MM_2004-2014_2-4.csv',
             'Lkb': './FLX_DE-Lkb_FLUXNET2015_FULLSET_2009-2013_1-4/FLX_DE-Lkb_FLUXNET2015_FULLSET_MM_2009-2013_1-4.csv',
             'Sor': './FLX_DK-Sor_FLUXNET2015_FULLSET_1996-2014_2-4/FLX_DK-Sor_FLUXNET2015_FULLSET_MM_1996-2014_2-4.csv',
             'Col': './FLX_IT-Col_FLUXNET2015_FULLSET_1996-2014_1-4/FLX_IT-Col_FLUXNET2015_FULLSET_MM_1996-2014_1-4.csv',
             'Ren': './FLX_IT-Ren_FLUXNET2015_FULLSET_1998-2013_1-4/FLX_IT-Ren_FLUXNET2015_FULLSET_MM_1998-2013_1-4.csv',
             'Fyo': './FLX_RU-Fyo_FLUXNET2015_FULLSET_1998-2014_2-4/FLX_RU-Fyo_FLUXNET2015_FULLSET_MM_1998-2014_2-4.csv'}


OBS_VARS = {"nee"  : ("NEE_VUT_REF", "kg m-2 month-1", "Net Ecosystem Exchange"), #fluxnet 
            "gpp"  : ("GPP_NT_VUT_REF", "kg m-2 month-1", "Gross Primary Productivity"),
            "reco" : ("RECO_NT_VUT_REF", "kg m-2 month-1", "Ecosystem Respiration"),
            "mle"  : ("LE_F_MDS", "W m-2", "Latent Heat Flux"), # COnvert to AET
            "tas"  : ("TA_F_MDS", "celcius", "air temperature"),
            "aet"  : ("AET", "kg m-2 month-1", "Actual Evapotranspiration")}


def get_conv_func(var):

    def id_func(input):
        return input

    if var == 'pr':
        return convert_pr
    elif var == 'ps':
        return convert_ps
    elif var == 'tas':
        return convert_ta
    elif var == 'vpd':
        return convert_vpd
    else:
        return id_func

def get_timestaps(site):
    data = pd.read_csv(OBS_SITES[site])["TIMESTAMP"].__array__()
    return f"{data[0]}01", f"{data[-1]}31"

def get_data(fpath, var=None):
    f = get_conv_func(var)
    return f(pd.read_csv(fpath)[FLUXNET_FULLSET_VARS[var][0]].__array__())

def get_ref_data(site, var=None):
    assert var in ['nee','gpp','reco', "tas"]
    if var == 'tas':
        return pd.read_csv(OBS_SITES[site])[OBS_VARS[var][0]].__array__()
    else:
        # need to convert C fluxes to make comparable        
        return pd.read_csv(OBS_SITES[site])[OBS_VARS[var][0]].__array__() * 0.0304368 # gm -2 d-1 => kg m-2 month-1

def calc_LHV(temp):
    f = np.vectorize(lambda P: 2.501 - (2.361 * 10e-3) * P)
    return f(temp)

def get_aet(site):
    mle = pd.read_csv(OBS_SITES[site])[OBS_VARS["mle"][0]].__array__() # W m-2
    tas = get_ref_data(site, "tas")
    mle *= 1e-6 # convert to MJ m-2 s-1
    return (mle / calc_LHV(tas)) * 2.62974e6 ## kg m-2 month-1

def create_arrs(var, mod_var=None):
    lat = []
    lon = []
    fdata = []
    names = []
    counter = 0
    for k, v in SITES_COORDINATES.items():
        # print(k, v)
        if k == 'SITE':
            continue
        counter += 1
        fpath = v[3]
        lat.append(v[0])
        lon.append(v[1])
        names.append(v[2])
        fdata.append(fpath)
    dt = get_data(fpath, var)

    out_data = np.zeros(shape=(counter, dt.size))

    for x in range(counter):
        out_data[x, :] = get_data(fdata[x], var)
    
           
    if mod_var is not None:
            out_data = mod_var(out_data)

    return out_data, np.concatenate(lat), np.concatenate(lon), np.array(names, dtype="<U7")

def create_gridlist(fname):
    with open(f"{fname}.txt", 'w', encoding="utf-8") as fh:
        pass
    with open(f"{fname}.txt", 'a', encoding="utf-8") as fh:
        for k, v in SITES_COORDINATES.items():
            # print(k, v)
            if k == 'SITE':
                continue
            lat = v[0][0]
            lon = v[1][0]
            fname = v[2]
            fh.write(f"{str(round(lon, 2))}\t{str(round(lat, 2))}\t{fname}\r\n")

def timeseries(fname = None,
              arr=None,
              var=None,
              unit=None,
              names = None,
              descr=None,
              time=None,
              la=None,
              lo=None,
              site_data=True,
              set="FULLSET",
              reference=None):

    """write fluxnet 2015 selected sites to nc"""


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
        dset.createDimension("station",size=arr.shape[0])
        dset.createDimension("time",size=arr.shape[1])

        # Data description
        dset.description = f"FLUXNET 2015 DATA - {set} {FLUXNET_FULLSET_VARS[var][0]}"
        dset.source = f'Forcing data for DVM - {descr}'
        dset.history= f'Created: {time_mod.ctime(time_mod.time())}'
        if reference is not None:
            dset.reference = reference
        dset.featureType = "timeSeries"

        # Create netCDF variables
        S = dset.createVariable("station", 'i4', ("station",), fill_value=999999)
        X  = dset.createVariable("lon", 'f4', ("station",), fill_value=1e+20)
        Y =  dset.createVariable("lat", 'f4', ("station",), fill_value=1e+20)
        SN = dset.createVariable("station_name", '<U7', ("station", ),fill_value= 'aaaaaaa')
        T  = dset.createVariable("time", 'i4', ("time",), fill_value=999999)

        D  = dset.createVariable(var, 'f4', ("station", "time"), fill_value=1e+20)

        S[...] = np.arange(arr.shape[0])
        T[...] = time_dim.__array__()
        T.units    = time_unit
        T.calendar = calendar
        T.standard_name = "time"
        T.axis = 'T'

        X[...] = lons
        X.units    = "degrees_east"
        X.long_name = 'station_longitude'
        X.standard_name = 'longitude'
        X.axis='X'

        Y[...] = lats
        Y.units    = "degrees_north"
        Y.long_name = 'station_latitude'
        Y.standard_name = 'latitude'
        Y.axis = 'Y'

        if names is None:
            nm = np.array(['aaaaaaa','aaaaaaa','aaaaaaa'], dtype='<U7')
        else:
            nm = names
        SN[...] = nm
        SN.long_name = "station name"
        SN.cf_role = "timeseries_id"


        D[...] = np.copy(arr[:,:], order="C")
        D.units = unit
        # D.fluxnet_name = FLUXNET_FULLSET_VARS[var][0]
        D.standard_name = FLUXNET_FULLSET_VARS[var][3]
        D.coordinates = u"lon lat"

        dset.close()

def cf_timeseries(fname = None,
                      arr=None,
                      var=None,
                      site=None,
                      unit=None,
                      descr=None,
                      time=None,
                      la=None,
                      lo=None,
                      set="FULLSET",
                      reference=None):

    """Create a CF compliant nc4 file for site(s) data TODO document
    Single Timeseries CF conventions (page 163)"""


    if fname is not None:
        dset = Dataset(os.path.join(Path('./'), fname + ".nc"),mode="w", format="NETCDF4")
    else:
        assert var is not None, "Need to set fname or var"
        dset = Dataset(os.path.join(Path('./'), OBS_VARS[var][0] + '.nc'), mode="w", format="NETCDF4")

    # Create temporal dimension
    time_dim = time['data']
    time_unit = time['units']
    calendar = time['calendar']

    lats = la
    lons = lo

    # Create netCDF dimensions
    dset.createDimension("time",size=arr.size)

    # Data description
    dset.description = f"FLUXNET 2015 DATA - {set} {OBS_VARS[var][2]}"
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

    X[...] = lons
    X.units    = "degrees_east"
    X.long_name = 'station_longitude'
    X.standard_name = 'longitude'

    Y[...] = lats
    Y.units    = "degrees_north"
    Y.long_name = 'station_latitude'
    Y.standard_name = 'latitude'

    SN[...] = SITES_COORDINATES[site][2]
    SN.long_name = "station name"
    SN.cf_role = "timeseries_id"


    D[...] = np.copy(arr[:], order="C")
    D.units = unit
    D.standard_name = var
    D.coordinates = u"time lat lon station_name"
    # D.coordinates = 'lon lat'

    dset.close()

def write_site_nc(VAR, mod=None):
    
    """write FLUXNET2015 FULLSET Variable to a netCDF4 file """
    start = "19890101"
    end= "20141231"
    idx = pd.date_range(start, end, freq='D')

    time_data = np.arange(idx.size, dtype='i4')

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
        arr, lat, lon, names = create_arrs(VAR)
        success = True
    except:
        if VAR == 'hurs':
            try:
                vpd, lat, lon, names = create_arrs('vpd')[0][0]
                # TODO error
                tair = create_arrs('tas')[0][0]
                arr = 100.0 * VPD2RH(tair, (vpd * (-1)) * 0.1)
                success = True
                assert False, "dont do hurs - error"
            except:
                success = False
        pass

    if success:
        la = lat
        lo = lon

        fname00 = f"{VAR}_FLUXNET2015"

        if mod is not None:
            ID, arr = mod(arr)
            fname00 = f"{VAR}_{ID}_FLUXNET2015"

        timeseries(fname=fname00, arr=arr, var=VAR, unit=units, names=names,
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

def write_ref_data(VAR, site):
    
    """write FLUXNET2015 FULLSET REFERENCE Variable to a netCDF4 file """
    
    assert VAR in ['nee', 'gpp', 'reco', "aet"]
    
    start, end = get_timestaps(site)
    idx = pd.date_range(start, end, freq='MS')
    
    # time_data = np.arange(idx.size, dtype='i4')

    day_init, hour_init = idx[0].isoformat().split("T")

    time_units = "days since %s %s" % (str(day_init), str(hour_init))

    calendar = 'proleptic_gregorian'
    
    time_data = cftime.date2num(idx.to_pydatetime(), units=time_units, calendar=calendar)
    
    descr = OBS_VARS[VAR][2]
    units = OBS_VARS[VAR][1]

    time_dict = {'data': time_data,
                'units': time_units,
                'calendar': calendar}

    if VAR == "aet":
        arr = get_aet(site)
    else:
        arr = get_ref_data(site, VAR)
    
    la = SITES_COORDINATES[site][0][0]
    lo = SITES_COORDINATES[site][1][0]
    
    fname00 = f"{VAR}_{site}_FLUXNET2015"

    cf_timeseries(fname=fname00, arr=arr, var=VAR, site=site,unit=units,
                    descr=descr, time=time_dict, la=la, lo=lo,
                    reference=FLUXNET_REF)

if __name__ == "__main__":
    pass
