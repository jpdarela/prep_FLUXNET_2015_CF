import numpy as np

Pa2kPa = 1e-3

#Converts precipitation from mm/day (FLUXNET) to kg m-2 s-1 (SI) :: assume 1 mm as 1 kg m-2
convert_pr = np.vectorize(lambda P: P * 1.5741e-5)

# VPD from hPa (FLUXNET) to kPa (GUESS) 
convert_vpd = np.vectorize(lambda P: (P * 0.1)) # * (-1))

# PS from kPa (FLUXNET) to Pa (GUESS)
convert_ps = np.vectorize(lambda P: P * 1e3)

convert_ta = np.vectorize(lambda P: P + 273.15)


def Vsat_slope(Tair:np.array,method=3) -> np.array:
# Translated to python from the bigleaf R package

#' Saturation Vapor Pressure (Esat) and Slope of the Esat Curve
#'    
#' @references Sonntag D. 1990: Important new values of the physical constants of 1986, vapor 
#'             pressure formulations based on the ITS-90 and psychrometric formulae. 
#'             Zeitschrift fuer Meteorologie 70, 340-344.
#'             
#'             World Meteorological Organization 2008: Guide to Meteorological Instruments
#'             and Methods of Observation (WMO-No.8). World Meteorological Organization,
#'             Geneva. 7th Edition.
#'             
#'             Alduchov, O. A. & Eskridge, R. E., 1996: Improved Magnus form approximation of 
#'             saturation vapor pressure. Journal of Applied Meteorology, 35, 601-609
#'             
#'             Allen, R.G., Pereira, L.S., Raes, D., Smith, M., 1998: Crop evapotranspiration -
#'             Guidelines for computing crop water requirements - FAO irrigation and drainage
#'             paper 56, FAO, Rome.

    """ Tair °C """
    
    methods = ("Sonntag_1990","Alduchov_1996","Allen_1998")  
    assert method <= 3 or method >= 1, "Methods:\n1 - Sonntag_1990\n2 - Alduchov_1996\n3 - Allen_1998"
    formula = methods[method - 1]
    print(f"Calculating RH from (TAIR, VPD) with Esat slope using the formula: {formula}")
    
    if formula == "Sonntag_1990":
        a = 611.2
        b = 17.62
        c = 243.12
    elif (formula == "Alduchov_1996"):
        a = 610.94
        b = 17.625
        c = 243.04
    elif (formula == "Allen_1998"):
        a = 610.8
        b = 17.27
        c = 237.3
  
  # saturation vapor pressure
    Esat = a * np.exp((b * Tair) / (c + Tair))
    Esat = Esat * Pa2kPa
  
  # slope of the saturation vapor pressure curve
    Delta = a * (np.exp((b * Tair)/(c + Tair)) * (b/(c + Tair) - (b * Tair)/(c + Tair)**2))
    Delta = Delta * Pa2kPa
  
    return Esat,Delta


def VPD2RH(Tair:np.array, VPD:np.array) -> np.array:
    """Estimate hurs from Tair (°C) and VPD (kPa)"""
    # Translated to python from the bigleaf R package
    esat =  Vsat_slope(Tair)[0]    
    return 1.0 - (VPD / esat)


