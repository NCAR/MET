#!/bin/sh

echo
echo "*** Running MODE on APCP using netCDF input for both forecast and observation ***"
../bin/mode \
   ../out/pcp_combine/sample_fcst_12L_2005080712V_12A.nc \
   ../out/pcp_combine/sample_obs_2005080712V_12A.nc \
   config/WrfModeConfig_APCP_12 \
   -outdir ../out/mode -v 2

echo
echo "*** Running MODE on APCP using a GRIB forecast and netCDF observation ***"
../bin/mode \
   ../data/sample_fcst/2005080700/wrfprs_ruc13_24.tm00_G212 \
   ../out/pcp_combine/sample_obs_2005080800V_24A.nc \
   config/WrfModeConfig_APCP_24 \
   -outdir ../out/mode -v 2

echo
echo "*** Running MODE on RH at 500mb using GRIB input for both forecast and observation ***"   
../bin/mode \
   ../data/sample_fcst/2005080700/wrfprs_ruc13_12.tm00_G212 \
   ../data/sample_fcst/2005080712/wrfprs_ruc13_00.tm00_G212 \
   config/WrfModeConfig_RH \
   -outdir ../out/mode -v 2
