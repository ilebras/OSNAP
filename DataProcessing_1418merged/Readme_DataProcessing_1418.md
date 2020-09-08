## This document briefly describes the scripts in this folder and outlines the steps taken to grid the first four years of Cape Farewell OSNAP data uniformly

Note: these scripts call my new auxiliary functions .py file -- firstfuncs_1618.py

All figures are output to figures_1418merged folder

Plots and outcomes are described in /OSNAP/data_processing_records/Merged_1418_dataprocessing/Data_Processing_Overview

## "FirstLook" scripts

* are the first look at the 2016-2018 microcat data -- take a peak at the higher frequencies and take a look at instrument placement
* Also, in FirstLook_TSdata_Grid_2018reco.py, I produce the daily averaged mcat files for 2018 recovery data

## "UnifyFormat" scripts

* Get every mooring/deployment into the same daily averaged format. Still instrument by instrument. The velocity data are ALL lowpass filtered (2nd order butterworth, 40hour cutoff). 
* Output from this stage is stored in netcdf files in Daily_netcdf folders within each 2016recovery/2018recovery folders
* Note that the filtering/dropping of bins/quick QC is done at this stage for velocity

## "Mooring_merged_corrections/checkout/reconstruct" scripts

* In these scripts I take a look at the daily microcat data and correct/remove as needed.
* Then they get saved as new files within Daily_netcdf folder with appropriate corr/recon addition.

## "GridVertically" scripts

* Here I get all the data onto a uniform 2m grid (note the switch to depth)
* For mcat data -- extend to surface using top 20m gradient in psal and ptmp -- then derive density. Extend bottomost measurement downward as needed.
* For vel data -- extend top measurement to surface, bottommost to deepest.
* The two occupations are merged at this point - output is stored in data/OSNAP_CFgridded_2014-2018

## "Combine_mcat_vel_xarrays.py"

* Load all mcat and velocity data, merge the two together, create a new netcdf file which includes all moorings, T,S and Vel

## "FinerGrid_AlongTopo_PlusShelf.py"

* Produces the final netcdf which includes:
  * extrapolation 10km onto the shelf based on horizontal gradient between CF1 and CF2, and geostrophy
  * Filling in of bottom triangles via interpolation along the bottom
  * All interpolation is linear.
  * Note: new horizontal grid includes 3 additional points between each pair of moorings (including the synthetic one at -10km).



