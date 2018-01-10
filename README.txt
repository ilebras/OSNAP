Log of hydrogen scripts and their functions

* aux_funcs.py: commonly loaded set of functions and classes to loaded

* InspectCorr_salinity_issues.py: script in which dictionary of CTD data were loaded and re-corrected, based on inversions and salinity drifts.
Note: temperature changed to potential temperature here! (python 2)

* InspectCorr_salinity_issues_15min.py: Implementation of corrections to 15 minute resolution

* Makenetcdf_salcorr.py

* Interp_CTD.py: CTD data loaded (from InspectCorr_salinity_issues where appropriate) and interpolated vertically. Tidbit data also optionally incorporated here (python 2)

* Interp_CTD_notid_pden.py: Update which corrects the implementation of and takes tidbits out of the interpolation for now.

* Inspect_saltmpinterp.py: Make some additional quick plots of newly interpolated data to see how it looks with and without tidbits, in TS space, and re: inversions.

* run_CTD.py: so far, a for loop to put all moorings through Interp_CTD

* Transport.py: so far, create xarray and get EGC system (freshwater) transport

* Makenetcdf_salcorr.py: overwrites salinity variable in netcdf files to the corrected value

* create_xarray.py: Gets saltmppden and vel on same grid and makes an xarray

* create_finergrid_xarrays.py: Uses output of create_xarray.py to make finer grid data directly usable for calculating transport and making plots

* Shipboard_CTD.py:  load Shipboard CTD data from several cruises and make xarray. Also a lot of auxiliary analysis in this script, including comparison with Nuka  re: surface extrapolate

* Shipboard_ADCP.py: load vessel mounted adcp data, plot and compare

* Shipboard_LADCP.py: load LADCP data and make xarray

* Shipboard_xport_LADCP_CTD.py: load xarrays and combine CTD and LADCP data

* Sections.py: somewhat of a legacy script, but still useful methods to draw upon in there.
