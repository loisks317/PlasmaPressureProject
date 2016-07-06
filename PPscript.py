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
sats=['B']
lsats=['b']
date1='20130601'
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
MLTLBins=[[[ [] for i in range(nL)] for j in range(nMLT)] for k in range(9)]
# 1 = < 1 up until 9 
#
# get the HOPE pressure
for iSat in range(len(sats)):
 dt1=dt0
 while dt1 < dt2:
  monthCur=dt1.month
  date1=datetime.datetime.strftime(dt1, '%Y%m%d')
  try:
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_MOM_'+sats[iSat])
    f='rbsp'+lsats[iSat]+'_rel03_ect-hope-MOM-L3_'+date1+'*.cdf'
    file=glob.glob(f)
    pyf=pycdf.CDF(file[0])
    #
    # get the RBSPICE data
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/RBSPICE_MOM_'+sats[iSat])
    fspice='rbsp-'+lsats[iSat]+'-rbspice_lev-3-PAP_TOFxEH_'+date1+'*.cdf'
    fileSpice=glob.glob(fspice)
    pyfspice=pycdf.CDF(fileSpice[0])
        
    #
    # get the pressures
    # first HOPE
    kelvin=11604.505
    boltzmann=1.38*1e-23
    hTime=pd.DatetimeIndex(pyf['Epoch_Ion'][...])
    density=np.array(pyf['Dens_p_30'][...])*1.0e6
    temperature=1/3.*np.array(pyf['Tpar_p_30'][...])+2/3.*np.array(pyf['Tperp_p_30'])
    HOPEpressure=(density*boltzmann*temperature)/1.0e-9 # to get it in nPa
    #
    # now put into pandas
    rt = pd.period_range(date1,periods=1440, freq='T').to_timestamp()
    rng=rt[::1]
    HPress=(pd.DataFrame(HOPEpressure, index=hTime, columns=['Pressure'])).resample('1min', how='median').reindex(index=rng)
    #
    # now RBSPICE
    RBTime=pd.DatetimeIndex(pyfspice['Epoch'][...])
    RBPressure=np.array(pyfspice['FPDU_ParaPressure'][...])/3. + 2*np.array(pyfspice['FPDU_PerpPressure'][...]/3.)
    RBPress= (pd.DataFrame(RBPressure, index=RBTime, columns=['Pressure'])).resample('1min', how='median').reindex(index=rng)
    TotalPressure=np.array(RBPress) + np.array(HPress)
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
                   MLTLBins[iKp][iMLT][iL]+=list(TotalPressure[matches])
    #
    # pickle by month
  except:
      print("Bad Data for " + date1)
  yrCur=dt1.year
  dt1=dt1+datetime.timedelta(days=1)
  if dt1.month !=monthCur:
     print(str(yrCur)+'_'+str(monthCur))
     os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaPressure/TotalPressure/')
     for iKp in range(9):
             pickle.dump(MLTLBins[iKp], open('Pressure_Kp='+str(iKp)+'_date='+str(yrCur)+'_'+str(monthCur)+'sat='+sats[iSat], 'wb'))
     monthCur=dt1.month
     MLTLBins=[[[ [] for i in range(nL)] for j in range(nMLT)] for k in range(9)]
