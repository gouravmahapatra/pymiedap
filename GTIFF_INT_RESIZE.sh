#!/bin/bash

obs_year="2011"
obs_type="MYD08_D3.A"
EXT=".tif"
directory=/Users/Ashwyn/Desktop/observation_database/

LC=$directory"MCD12Q1.006_LC_Type1_doy2011001_aid0001.tif"
LCnew=$directory"MCD12Q1.006_LC_Type1_doy2011_SAMP.tif"

gdal_translate -of GTiff -outsize 720 360 -r mode $LC $LCnew

for filename in $directory"$obs_type$obs_year"*Cloud_Fraction_Mean.tif; do
    date=$(echo $filename| cut -d'.' -f 2)
    prefix="A$obs_year"
    day=${date#"$prefix"}

    CF=$directory"$obs_type$obs_year$day*Cloud_Fraction_Mean.tif"
    CTP=$directory"$obs_type$obs_year$day*Cloud_Top_Pressure_Mean.tif"
    COT=$directory"$obs_type$obs_year$day*Cloud_Optical_Thickness_Liquid_Log_Mean.tif"
    CER=$directory"$obs_type$obs_year$day*Cloud_Effective_Radius_Liquid_Mean.tif"
    OW1=$directory"$obs_type$obs_year$day.00-Ocean_Surface_Wind_Speed_10m.tif"
    OW2=$directory"$obs_type$obs_year$day.06-Ocean_Surface_Wind_Speed_10m.tif"
    OW3=$directory"$obs_type$obs_year$day.12-Ocean_Surface_Wind_Speed_10m.tif"
    OW4=$directory"$obs_type$obs_year$day.18-Ocean_Surface_Wind_Speed_10m.tif"

    CFnew=$directory$obs_type$obs_year$day"_Cloud_Fraction_Mean_INT_SAMP.tif"
    CTPnew=$directory$obs_type$obs_year$day"_Cloud_Top_Pressure_Mean_INT_SAMP.tif"
    COTnew=$directory$obs_type$obs_year$day"_Cloud_Optical_Thickness_Liquid_Log_Mean_INT_SAMP.tif"
    CERnew=$directory$obs_type$obs_year$day"_Cloud_Effective_Radius_Liquid_Mean_INT_SAMP.tif"
    OWnew1=$directory$obs_type$obs_year$day".00-Ocean_Surface_Wind_Speed_10m_INT_SAMP.tif"
    OWnew2=$directory$obs_type$obs_year$day".06-Ocean_Surface_Wind_Speed_10m_INT_SAMP.tif"
    OWnew3=$directory$obs_type$obs_year$day".12-Ocean_Surface_Wind_Speed_10m_INT_SAMP.tif"
    OWnew4=$directory$obs_type$obs_year$day".18-Ocean_Surface_Wind_Speed_10m_INT_SAMP.tif"

    gdal_calc.py -A $CF --overwrite --NoDataValue=0 --outfile=CFmask.tif --calc="1*(A<0)"
    gdal_calc.py -A $CTP --overwrite --NoDataValue=0 --outfile=CTPmask.tif --calc="1*(A<0)"  #1-> int, 0-> no int
    gdal_calc.py -A $COT --overwrite --NoDataValue=0 --outfile=COTmask.tif --calc="1*(A<0)"
    gdal_calc.py -A $CER --overwrite --NoDataValue=0 --outfile=CERmask.tif --calc="1*(A<0)"


    gdal_calc.py -A CFmask.tif --overwrite --NoDataValue=1 --outfile=CFmask.tif --calc="1*(A==0)"
    gdal_calc.py -A CTPmask.tif --overwrite --NoDataValue=1 --outfile=CTPmask.tif --calc="1*(A==0)"  #0-> int, 1-> no int
    gdal_calc.py -A COTmask.tif --overwrite --NoDataValue=1 --outfile=COTmask.tif --calc="1*(A==0)"
    gdal_calc.py -A CERmask.tif --overwrite --NoDataValue=1 --outfile=CERmask.tif --calc="1*(A==0)"

    gdal_fillnodata.py -md 10 $CF -mask CFmask.tif CFnew.tif -of GTiff
    gdal_fillnodata.py -md 10 $CTP -mask CTPmask.tif CTPnew.tif -of GTiff
    gdal_fillnodata.py -md 10 $COT -mask COTmask.tif COTnew.tif -of GTiff
    gdal_fillnodata.py -md 10 $CER -mask CERmask.tif CERnew.tif -of GTiff

    gdal_calc.py -A CTPnew.tif -B $CF --overwrite --NoDataValue=0 --outfile=CTPnew.tif --calc="-9999*(B==0)+A*(B<0)+A*(B>0)"
    gdal_calc.py -A COTnew.tif -B $CF --overwrite --NoDataValue=0 --outfile=COTnew.tif --calc="-9999*(B==0)+A*(B<0)+A*(B>0)"
    gdal_calc.py -A CERnew.tif -B $CF --overwrite --NoDataValue=0 --outfile=CERnew.tif --calc="-9999*(B==0)+A*(B<0)+A*(B>0)"

    echo $OW1
    gdal_translate -of GTiff -outsize 720 360 -r nearest CFnew.tif $CFnew
    gdal_translate -of GTiff -outsize 720 360 -r nearest CTPnew.tif $CTPnew
    gdal_translate -of GTiff -outsize 720 360 -r nearest COTnew.tif $COTnew
    gdal_translate -of GTiff -outsize 720 360 -r nearest CERnew.tif $CERnew
    gdal_translate -of GTiff -outsize 720 360 -r nearest $OW1 $OWnew1
    gdal_translate -of GTiff -outsize 720 360 -r nearest $OW2 $OWnew2
    gdal_translate -of GTiff -outsize 720 360 -r nearest $OW3 $OWnew3
    gdal_translate -of GTiff -outsize 720 360 -r nearest $OW4 $OWnew4

done
rm CERmask.tif
rm CERnew.tif
rm CFmask.tif
rm CFnew.tif
rm COTmask.tif
rm COTnew.tif
rm CTPmask.tif
rm CTPnew.tif
