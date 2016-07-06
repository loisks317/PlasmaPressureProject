# loadPP.py
#
# load in the Plasma Pressure and plot by Kp
#
# LKS SANSA January 2016
#
# imports
import numpy as np
import glob
import os
import pickle
import datetime
import matplotlib.pyplot as plt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
import matplotlib.collections as collections
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter 
import matplotlib.dates as dt
from dateutil.relativedelta import relativedelta
from matplotlib.colors import LogNorm
from numpy import ma
os.chdir('/Users/loisks/Desktop/Functions/')
import colormaps as cmaps
os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaPressure/')
plt.register_cmap(name='viridis', cmap=cmaps.viridis)

def dual_half_circle(center, radius, angle=90, ax=None, colors=('w','k'),
                     **kwargs):
    from matplotlib.patches import Wedge
    """
    Add two half circles to the axes *ax* (or the current axes) with the 
    specified facecolors *colors* rotated at *angle* (in degrees).
    """
    if ax is None:
        ax = plt.gca()
    theta1, theta2 = angle, angle + 180
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[0],transform=ax.transData._b, **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], transform=ax.transData._b,**kwargs)
    for wedge in [w1, w2]:
        ax.add_artist(wedge)
    return [w1, w2]

def fancy_plot(data, mltbins, Lbins, vmin, vmax, colorbar_Label, species, activity, plot_label, directory, scale,cmap) :
    # data = your data, for this plot it needs to be a 2-D array
    # mltbins = # of mltbins
    # Lbins = # of Lbins
    # vmin = min of data, currently set to log
    # vmax = max of data, currently set to log
    # colorbar_label = label for colorbar
    # species = H, He, or O? 
    # activity = high or low Kp
    # plot_label = what you want to name the files
    # directory = name of directory to store the files
    # scale = set to 'log' or 'linear'
    mlt=(np.linspace(0, 24, mltbins))*(15*np.pi/180.)
    L_shell=np.linspace(1.5,7.25 , Lbins)
    # intitialize the figure
    fig=plt.figure()
    ax=fig.add_subplot(111, polar=True)
    datah_m=ma.masked_invalid(np.array(data).transpose())
    X,Y=np.meshgrid(mlt, L_shell)
    #ax.set_ylim(0, 6.25)
    ax.set_ylim(0,7.25)
    #plt.title(title, fontsize=30, y=1.15)
    dual_half_circle((0,0), 1.0, angle=90, ax=None, colors=('w','k'))
    cbaxes = fig.add_axes([0.75, 0.15, 0.03, 0.75])
    if scale=='log':
        col=ax.pcolormesh( X, Y, np.log10(datah_m), cmap=cmap, vmin=vmin, vmax=vmax )
        cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1))
    else:
        col=ax.pcolormesh( X, Y,datah_m, cmap='jet', vmin=vmin, vmax=vmax)
        cb = plt.colorbar(col, cax = cbaxes,ticks=range(vmin, vmax+1))
    plt.subplots_adjust(right=0.73, top=0.8, bottom=0.15)
    xT2=['', '2', '', '4', '', '6', '']
    #xT2=['', '1', '', '2', '', '3']
    xL2=['00','', '06','', '12','', '18']
    cb.set_label(colorbar_Label, fontsize=30)
    ax.tick_params(axis='both', which='major', labelsize=25)    
    cb.ax.tick_params(labelsize=35) 
    ax.set_yticklabels(xT2, fontsize=30)
    ax.set_xticklabels(xL2, fontsize=30)
    ax.grid(True)
    plt.draw()
    subdir_name=directory
    if not os.path.exists(subdir_name):
        os.umask(0) # unmask if necessary
        os.makedirs(subdir_name) 
    os.chdir(subdir_name)#
    fig.set_size_inches(13,9)
    plt.savefig(species+'_'+activity+'_'+plot_label+'.png')
    plt.close(fig)
    os.chdir('..')

# satellites
sats=['A','B']
lsats=['a','b']
date1='20130201'
date2='20150401'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')

# bin sizes
nMLT=mltbins=48
nL=nLbins=24
MLTbins=np.linspace(0.25, 23.75, nMLT)
MLTarr=[ [] for i in range(nMLT)]
Larr=[ [] for i in range(nL)]
Lbins=np.linspace(1.5, 7.25, nL)
MLTLBins=[[[ [] for i in list(range(nL))] for j in list(range(nMLT))] for k in range(9)]
#
def diff_month(d1, d2):
        return (d1.year - d2.year)*12 + d1.month - d2.month
    #

total_month=diff_month(dt2, dt0)
for iKp in range(9):
 dt1=dt0 
 cfilee=[[ [] for x in list(range(nLbins))] for x in list(range(mltbins))]
 meanfe=[[ [] for x in list(range(nLbins))] for x in list(range(mltbins))]
 for imonth in range(total_month):
    cur_date=str(dt1.year)+'_'+str(dt1.month)
    dt1=dt1+relativedelta(months=1)
    os.chdir('TotalPressure')
    file1=glob.glob('Pressure_Kp='+str(iKp)+'_date='+cur_date+'sat=A')[0]
    file2=glob.glob('Pressure_Kp='+str(iKp)+'_date='+cur_date+'sat=B')[0]
    f=pickle.load(open(file1,'rb'))
    f2=pickle.load(open(file2,'rb'))
    os.chdir('..')
    for imlt in range(nMLT):
        for iL in range(nL):
    # combine A and B
                cfilee[imlt][iL]=list(cfilee[imlt][iL])+list(f[imlt][iL])+list(f2[imlt][iL])
 #
 # need to screen for bad pressure
 for imlt in range(nMLT):
    for iL in range(nL):
          cfilee[imlt][iL]=np.array(cfilee[imlt][iL])
          bad=np.where(cfilee[imlt][iL]>1000)[0]
          cfilee[imlt][iL][bad]=np.nan
          meanfe[imlt][iL]=np.nanmean(cfilee[imlt][iL])

#
 fancy_plot(meanfe, nMLT, nL, -1, 1, 'Pressure [nPa]', 'H$^{+}$', 'Kp='+str(iKp), 'TotalPressure', 'TotalPressurePlots', 'log','magma')
