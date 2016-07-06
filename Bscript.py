# PPscript.py
#
# Call different modules to calculate total Plasma Pressure
# in the Van Allen Probes. Interpolate to 1 minute intervals
# Combine HOPE and RB-SPICE pressure moments
#
# LKS, 2016 SANSA
#
#imports 
import os
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
import datetime
import glob
from matplotlib import dates as dates
from datetime import timedelta
from dateutil import rrule
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from spacepy import pycdf
import glob
from spacepy import datamodel as dm
from matplotlib.ticker import LogLocator, LogFormatter, MultipleLocator, FormatStrFormatter, FixedLocator, FixedFormatter
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
import pandas as pd
import pickle
#
# satellites
sats=['A','B']
lsats=['A','b']
date1='20130201'
date2='20150501'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
#
# load in the Kp
KpArr=pickle.load(open('KpFeb2013_Apr2015.p', 'rb'))
#
# bin sizes
nMLT=48
nL=24
MLTbins=np.linspace(0.25, 23.75, nMLT)
MLTarr=[ [] for i in range(nMLT)]
Larr=[ [] for i in range(nL)]
Lbins=np.linspace(1.5, 7.25, nL)
BxAll=[[[ [] for i in range(nL)] for j in range(nMLT)] for k in range(9)]
ByAll=[[[ [] for i in range(nL)] for j in range(nMLT)] for k in range(9)]
BzAll=[[[ [] for i in range(nL)] for j in range(nMLT)] for k in range(9)]
BAll=[[[ [] for i in range(nL)] for j in range(nMLT)] for k in range(9)]
# 1 = < 1 up until 9 
#
# get the HOPE pressure
for iSat in range(len(sats)):
 dt1=dt0
 while dt1 < dt2:
  monthCur=dt1.month
  date1=datetime.datetime.strftime(dt1, '%Y%m%d')
  try:
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMFISIS_GEO_'+sats[iSat])
    f='rbsp-'+lsats[iSat]+'_magnetometer_1sec-geo_emfisis-L3_'+date1+'*.cdf'
    file=glob.glob(f)
    pyfe=pycdf.CDF(file[0])
    #
    time_emf=pyfe['Epoch'][...]
    pdETime=pd.DatetimeIndex(time_emf)
    coordinates=np.swapaxes(pyfe['coordinates'][...],1,0)
    emf_x=np.swapaxes(pyfe['Mag'][...], 1, 0)[0]
    emf_y=np.swapaxes(pyfe['Mag'][...], 1, 0)[1]
    emf_z=np.swapaxes(pyfe['Mag'][...], 1, 0)[2]
# subtract out IGRF
    os.chdir('/Users/loisks/Documents/Functions/')
    import IGRF
    import IGRF_support as RBI
    spherical_coords=RBI.geo_to_lat(coordinates[0], coordinates[1], coordinates[2])
    BIx=np.zeros(len(emf_x)); BIy=np.zeros(len(emf_x)); BIz=np.zeros(len(emf_x))
    for a in range(len(emf_x)):
         igrf_res=IGRF.IGRF_fortran(spherical_coords[2][a],spherical_coords[1][a] ,spherical_coords[0][a],int(yrmo[0:4]))
         transform=IGRF.convert_IGRF(igrf_res[2], igrf_res[0],igrf_res[1], spherical_coords[1][a], spherical_coords[2][a])	      	
         BIx[a]=transform[0]
         BIy[a]=transform[1]
         BIz[a]=transform[2]
    cBx=np.array(emf_x)-np.array(BIx)
    cBy=np.array(emf_y)-np.array(BIy)
    cBz=np.array(emf_z)-np.array(BIz)
    #
    # pandas these into arrays
    # may need to change indexing because of IGRF
    rt = pd.period_range(date1,periods=1440, freq='T').to_timestamp()
    rng=rt[::1]
    GEOBx=(pd.DataFrame(cBx, index=pdETime, columns=['Bx'])).resample('1min', how='median').reindex(index=rng)
    GEOBy=(pd.DataFrame(cBy, index=pdETime, columns=['By'])).resample('1min', how='median').reindex(index=rng)
    GEOBz=(pd.DataFrame(cBz, index=pdETime, columns=['Bz'])).resample('1min', how='median').reindex(index=rng)
    #
    # interpolate the coordinates too
    GEOX=(pd.DataFrame(coordinates[0], index=pdETime, columns=['x'])).resample('1min', how='median').reindex(index=rng)
    GEOY=(pd.DataFrame(coordinates[1], index=pdETime, columns=['y'])).resample('1min', how='median').reindex(index=rng)
    GEOZ=(pd.DataFrame(coordinates[2], index=pdETime, columns=['z'])).resample('1min', how='median').reindex(index=rng)
    #
    #
    desired_coords=['GSE']
    for ic in range(len(desired_coords)):
        Bx=np.zeros(1440); By=np.zeros(1440); Bz=np.zeros(1440)
        XGSE=np.zeros(1440); YGSE=np.zeros(1440); ZGSE=np.zeros(1440)
        for igp in range(1440):
         # next step here is to do this all in pandas         
               ctime=rng[ic]
               year=ctime.year
               day=dt1.timetuple()[7] # day number in year
               hour=ctime.hour
               minute=ctime.minute
               second=ctime.second
# now apply each coordinate system

               cvals=coord.Coords([[GEOX[igp],GEOY[igp],GEOZ[igp]]], 'GEO', 'car')
               bcvals=coord.Coords([[GEOBx[igp],GEOBy[igp],GEOBz[igp]]], 'GEO', 'car')
               cvals.ticks=Ticktock([ctime], 'UTC')
               bcvals.ticks=Ticktock([ctime], 'UTC')
               new_coord=cvals.convert(desired_coords[ic],'car')
               newB_coord=bcvals.convert(desired_coords[ic], 'car')
               XGSE[igp]=new_coord.x; YGSE[igp]=new_coord.y; ZGSE[igp]=new_coord.z
               Bx[igp]=newB_coord.x; By[igp]=newB_coord.y; Bz[igp]=newB_coord.z

            
    #
    # now we have the interpolated pressures... get the Kp, L, and MLT
    # from HOPE Ephemeris
    #
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMPHEMERIS_'+sats[iSat])
    f2=glob.glob('rbsp'+lsats[iSat]+'_def_MagEphem_OP77Q_'+date1+'*.txt')[0]
    pyf2=dm.readJSONheadedASCII(f2)
    # get params
    LShell=np.array(np.nanmean(np.array(pyf2['L']), axis=1)) # weird multiple rows with L
    MLT=np.array(pyf2['CDMAG_MLT'])
    Kp=np.array(KpArr[date1]) # in minutes 
    #
    # now we sort and compare
    # 0.5 MLT, 0.25 L, and Kp =1 bins up to kp 9   
    for iKp in range(9):
        for iMLT in range(nMLT):
            for iL in range(nL):
                 
                   temp=np.where((MLT >= MLTbins[iMLT]-0.25) & (MLT < MLTbins[iMLT]+0.25))[0]
                   temp2=np.where((LShell >= Lbins[iL]-.125) & (LShell < Lbins[iL]+0.125))[0]
                   temp3=np.where(( Kp >= iKp) & ( Kp < iKp+1))[0]
                   # get the set
                   matches=list(set(temp) & set(temp2) & set(temp3))
                   BxAll[iKp][iMLT][iL]+=list(Bx[matches])
                   ByAll[iKp][iMLT][iL]+=list(By[matches])
                   BzAll[iKp][iMLT][iL]+=list(Bz[matches])
                   BAll[iKp][iMLT][iL]+=list(np.sqrt(np.array(Bx[matches]**2) + np.array(By[matches]**2) + np.array(Bz[matches]**2))) 
    #
    # pickle by month
  except:
      print("Bad Data for " + date1)
  yrCur=dt1.year
  dt1=dt1+datetime.timedelta(days=1)
  if dt1.month !=monthCur:
     print(str(yrCur)+'_'+str(monthCur))
     os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaPressure/BData/')
     for iKp in range(9):
             pickle.dump(BxAll[iKp], open('Bx_Kp='+str(iKp)+'_date='+str(yrCur)+'_'+str(monthCur)+'sat='+sats[iSat], 'wb'))
             pickle.dump(ByAll[iKp], open('By_Kp='+str(iKp)+'_date='+str(yrCur)+'_'+str(monthCur)+'sat='+sats[iSat], 'wb'))
             pickle.dump(BzAll[iKp], open('Bz_Kp='+str(iKp)+'_date='+str(yrCur)+'_'+str(monthCur)+'sat='+sats[iSat], 'wb'))
             pickle.dump(BAll[iKp], open('B_Kp='+str(iKp)+'_date='+str(yrCur)+'_'+str(monthCur)+'sat='+sats[iSat], 'wb'))
     monthCur=dt1.month
     BxAll=[[[ [] for i in range(nL)] for j in range(nMLT)] for k in range(9)]
     ByAll=[[[ [] for i in range(nL)] for j in range(nMLT)] for k in range(9)]
     BzAll=[[[ [] for i in range(nL)] for j in range(nMLT)] for k in range(9)]
     BAll=[[[ [] for i in range(nL)] for j in range(nMLT)] for k in range(9)]
