from pylab import *
import gsw

# %% Code to process Aanderaa oxygen optode data integrated into the CTD
# % Prepared by H. Palevsky, based on code used for OOI Irminger 6, Aug. 2019
# Translated into python on AR45 OSNAP cruise by Isabela Le Bras on June 26, 2020


def aaoptode_sternvolmer(foil_coef, optode_phase, tempc, salin, press):
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %% AAOPTODE_STERNVOLMER: For Aanderaa optode 4XXX series, to interpret
    # % calphase measurements based on the Stern Volmer equation. To be used in
    # % combination with aaoptode_salpresscorr, which handles salinity and
    # % pressure corrections.
    #
    # %% INPUTS
    # % foil_coef:    Struct of calibration coefficients
    # % optode_phase: vector of optode phase measurements (calphase or dphase)
    # % tempc:        temperature in deg C
    # % salin:        Salinity
    # % press:        pressure (db) for pressure correction
    #
    # %% OUTPUTS
    # % optode_uM: Oxygen concentration in uM
    # % optode_umolkg: Oxygen concentration in umol/kg
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % H. Palevsky, 8/6/2019, based on sg_aaoptode.m


    # % Uchida et al. 2008 Stern-Volmer based calbration mode
    C = foil_coef;
    Ksv = C[0] + C[1]*tempc+C[2]*tempc**2;
    P0 = C[3] + C[4]*tempc;
    PC = C[5] + C[6]*optode_phase;
    optode_uM = ((P0/PC)-1)/Ksv;
    # % note - this uses SR, reference salinity
    SA = 35.16504/35*salin;
    CT = gsw.CT_from_t(SA,tempc,press);
    optode_umolkg = 1000*optode_uM/(1000+ gsw.sigma0(SA,CT));
    return optode_uM, optode_umolkg

def aaoptode_salpresscorr(O2raw,temp,sal,press,S0):
    # % AAOPTODE_SALPRESSCORR:  Calculates oxygen concentration in uM from raw O2
    # % concentration, correcting for temperature and salinity. Based on Salinity
    # % and Depth compensation equations in the Aanderaa operating manual (also
    # % see function sg_aaoptode.m)
    # %
    # % INPUTS:
    # % O2raw:       O2 concentration output from optode (in uM)
    # % temp:        temperature in deg C
    # % sal:         salinity (measured)
    # % press:       pressure (db) for pressure correction
    # % S0:          internal salinity setting for optode (default is zero)
    # %
    # % OUTPUT:
    # % O2corr:       Oxygen concentration in uM
    # % O2sat:        Oxygen percent saturation (%) (100 = saturated)
    #
    #
    temps = log((298.15-temp)/(273.15+temp)); #%scaled temperature
    SB = [-6.24097E-3, #%salinity correction coefficients (B0, B1, B2, and B3)
            -6.93498E-3,
            -6.90358E-3,
            -4.29155E-3];
    SC = -3.11680E-7; #%salinity correction coefficient C0
    D = 0.032; #%depth correction coefficient

    O2corr = O2raw*exp((sal - S0)*(SB[0]+SB[1]*temps+SB[2]*temps**2+SB[3]*temps**3)+ SC*(sal**2 - S0**2))*(1+press*D/1000);
    return O2corr

def aaoptode_ctdconvert(salin,tempc,press,volts):
# %% INPUT
# % volts: raw optode output from CTD in volts
# % tempc: temperature (deg C) from CTD
# % salin: salinity from CTD
# % press: pressure from CTD
#
# %THESE COEFFICIENTS ARE DIFFERENT FOR EACH OPTODE - they are found on the
# %calibration certificate sheet
# %You will need to update for specific sensor used - this is for SN277 (AR45 OSNAP cruise)
# %You should also check the salinity setting for the sensor - if it is
# %something other than zero, the final line of this script will need to be
# %modified

    # foilcoeff = array([2.82567E-3, 1.20716E-4, 2.4593E-6, 2.30757E2, -3.09502E-1, -5.60627E1, 4.5615E0]) #from Hilary
    # conccoeff = array([-1.28596, 1.039998]) #from Hilary
    foilcoeff =array([2.798512E-03, 1.179460E-04, 2.512907E-06, 2.262806E+02, -3.570254E-01, -6.104725E+01, 4.558537E+00]) #from Roo
    conccoeff  = array([0,   1.16]) #from Roo
    # % Pressure compensation coefficient should be 0.032 for all 4XXX sensors
    D = 0.032;
    ### Note: I took out the optional test of the 2-point calibration since I don't have that info.

    # %% Coefficients for processing data from Aanderaa optode
    # %Coefficients for converting from voltage to calphase
    A = 10;
    B = 12; #%note that this should be the same for all optodes
    calphase = B*volts + A;

    # %% Read data in from  CTD casts
    [optode_uM_preconn, na] = aaoptode_sternvolmer(foilcoeff, calphase, tempc, salin, press);
    optode_uM = conccoeff[0] + conccoeff[1]*optode_uM_preconn;
    O2corr = aaoptode_salpresscorr(optode_uM, tempc, salin, press, 0);
    return O2corr


##### Example data from Hilary:
tempc=array([5.9207,2.6164])
salin=array([ 34.7238,35.1903])
press=array([51,2301])
volts=array([1.956,2.2699])
O2corr_ex=array([373.605,337.0803])

O2corr=aaoptode_ctdconvert(salin,tempc,press,array([NaN,NaN]))

#verify that it matches what Hilary provided:
O2corr
