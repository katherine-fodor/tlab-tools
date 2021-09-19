#!/usr/bin/env python
# coding: utf-8

# In[2]:


#!/usr/bin/python3

from ReadStats import Statistics, Pdfs 
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
from matplotlib import rc

rc('text', usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font', family='serif')
rc('font', size=24)
rc('axes', linewidth=1.5)
rc('axes', labelsize=24)
rc('lines', linewidth=2)

opath = '/Volumes/Seagate/SCRATCH/plots/3D/Re025/Rapids/'

path_S0    = '/Volumes/Seagate/SCRATCH/qCBL_3D/Re025/2496x512x2496-S00/'
path_S05   = '/Volumes/Seagate/SCRATCH/qCBL_3D/Re025/1280x512x1280-S05/'
path_S10   = '/Volumes/Seagate/SCRATCH/qCBL_3D/Re025/1280x512x1280-S10/'
path_S15   = '/Volumes/Seagate/SCRATCH/qCBL_3D/Re025/1536x576x1536-S15/'
path_S20_1 = '/Volumes/Seagate/SCRATCH/qCBL_3D/Re025/S20-1536x576x1536/'
path_S20_2 = '/Volumes/Seagate/SCRATCH/qCBL_3D/Re025/S20-1536x576x2304/'
path_S25   = '/Volumes/Seagate/SCRATCH/qCBL_3D/Re025/2560x896x2560-S25/'


# In[4]:


# Constants

nu = 1./15000.
B0 = 0.005
N = np.sqrt(3)
L0 = (B0/N**3)**0.5
ceps = 0.1


# In[3]:


def runningmean(timeseries,window):
    """
    Calculate a running mean of a time series with a specified window size.
    Window size is number of entries on either side of the entry being averaged.
    Leaves out entries at the beginning and end of the time series such that the
    window size is always the same, but the resulting time series is shorter than
    the original.
    """
    
    nt = len(timeseries)
    outseries = np.zeros(nt-(window*2))
    for n in range(window,nt-window):
        outseries[n-window] = np.mean(timeseries[n-window:n+window+1])
    return outseries


# In[4]:


#Stats
# Data is at every one z_enc/L_0 from z_enc/L_0 = 15 to 30.

S0    = Statistics(path_S0+'stats/pdftimes/avg9500-53000.nc')
S05   = Statistics(path_S05+'stats/pdftimes/avg12000-67000.nc')
S10   = Statistics(path_S10+'stats/pdftimes/avg13000-84000.nc')
S15   = Statistics(path_S15+'stats/pdftimes/avg15000-92000.nc')
S20_1 = Statistics(path_S20_1+'stats/pdftimes/avg17000-35000.nc')
S20_2 = Statistics(path_S20_2+'stats/pdftimes/avg39000-91000.nc')
S25   = Statistics(path_S25+'stats/pdftimes/avg28000-128000.nc')

z_enc_S20 = np.concatenate((S20_1.z_enc,S20_2.z_enc))
z_enc_arg_S20 = np.concatenate((S20_1.z_enc_arg,S20_2.z_enc_arg))

z_enc_runningmean_S0 = runningmean(S0.z_enc,1)
z_enc_runningmean_S05 = runningmean(S05.z_enc,1)
z_enc_runningmean_S10 = runningmean(S10.z_enc,1)
z_enc_runningmean_S15 = runningmean(S15.z_enc,1)
z_enc_runningmean_S20 = runningmean(z_enc_S20,1)
z_enc_runningmean_S25 = runningmean(S25.z_enc,1)


# In[5]:


# Pdfs
# Data is at every one z_enc/L_0 from z_enc/L_0 = 15 to 30.

vortlist_S0 = [path_S0+'stats/pdfs/pdf9500.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf11000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf13000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf14500.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf16000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf18000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf21000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf23000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf26000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf29000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf32000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf36000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf40000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf44000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf48000.LnEnstrophyW_iW_i',
               path_S0+'stats/pdfs/pdf53000.LnEnstrophyW_iW_i']

vortlist_S05 = [path_S05+'stats/pdfs/pdf12000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf14000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf17000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf19000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf22000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf25000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf28000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf31000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf34000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf38000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf42000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf47000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf52000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf57000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf62000.LnEnstrophyW_iW_i',
                path_S05+'stats/pdfs/pdf67000.LnEnstrophyW_iW_i']

vortlist_S10 = [path_S10+'stats/pdfs/pdf13000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf16000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf18000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf21000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf24000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf28000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf32000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf37000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf41000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf46000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf52000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf58000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf64000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf70000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf77000.LnEnstrophyW_iW_i',
                path_S10+'stats/pdfs/pdf84000.LnEnstrophyW_iW_i']

vortlist_S15 = [path_S15+'stats/pdfs/pdf15000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf18000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf21000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf24000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf27000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf31000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf35000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf40000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf45000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf50000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf56000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf63000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf69000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf76000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf84000.LnEnstrophyW_iW_i',
                path_S15+'stats/pdfs/pdf92000.LnEnstrophyW_iW_i']

vortlist_S20_1 = [path_S20_1+'stats/pdfs/pdf17000.LnEnstrophyW_iW_i',
                  path_S20_1+'stats/pdfs/pdf19000.LnEnstrophyW_iW_i',
                  path_S20_1+'stats/pdfs/pdf22000.LnEnstrophyW_iW_i',
                  path_S20_1+'stats/pdfs/pdf25000.LnEnstrophyW_iW_i',
                  path_S20_1+'stats/pdfs/pdf28000.LnEnstrophyW_iW_i',
                  path_S20_1+'stats/pdfs/pdf31000.LnEnstrophyW_iW_i',
                  path_S20_1+'stats/pdfs/pdf35000.LnEnstrophyW_iW_i']
vortlist_S20_2 = [path_S20_2+'stats/pdfs/pdf39000.LnEnstrophyW_iW_i',
                  path_S20_2+'stats/pdfs/pdf44000.LnEnstrophyW_iW_i',
                  path_S20_2+'stats/pdfs/pdf50000.LnEnstrophyW_iW_i',
                  path_S20_2+'stats/pdfs/pdf63000.LnEnstrophyW_iW_i',
                  path_S20_2+'stats/pdfs/pdf64000.LnEnstrophyW_iW_i',
                  path_S20_2+'stats/pdfs/pdf67000.LnEnstrophyW_iW_i',
                  path_S20_2+'stats/pdfs/pdf75000.LnEnstrophyW_iW_i',
                  path_S20_2+'stats/pdfs/pdf84000.LnEnstrophyW_iW_i',
                  path_S20_2+'stats/pdfs/pdf91000.LnEnstrophyW_iW_i']

vortlist_S25 = [path_S25+'stats/pdfs/pdf28000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf32000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf36000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf41000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf46000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf51000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf57000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf63000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf70000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf77000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf85000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf93000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf101000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf110000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf119000.LnEnstrophyW_iW_i',
                path_S25+'stats/pdfs/pdf128000.LnEnstrophyW_iW_i']


vortpdf_S0    = Pdfs(vortlist_S0,path_S0+'y.dat')
vortpdf_S05   = Pdfs(vortlist_S05,path_S05+'y.dat')
vortpdf_S10   = Pdfs(vortlist_S10,path_S10+'y.dat')
vortpdf_S15   = Pdfs(vortlist_S15,path_S15+'y.dat')
vortpdf_S20_1 = Pdfs(vortlist_S20_1,path_S20_1+'y.dat')
vortpdf_S20_2 = Pdfs(vortlist_S20_2,path_S20_2+'y.dat')
vortpdf_S25   = Pdfs(vortlist_S25,path_S25+'y.dat')


# In[6]:


# Create grid on which to interpolate pdfs
# Doesn't really matter which one is used, as long as it's the same in PDFs.ipynb

S0_vortpdf_interp_data    = Pdfs([path_S0+'stats/pdfs/pdf48000.LnEnstrophyW_iW_i'],path_S0+'y.dat')
S05_vortpdf_interp_data   = Pdfs([path_S05+'stats/pdfs/pdf57000.LnEnstrophyW_iW_i'],path_S05+'y.dat')
S10_vortpdf_interp_data   = Pdfs([path_S10+'stats/pdfs/pdf84000.LnEnstrophyW_iW_i'],path_S10+'y.dat')
S15_vortpdf_interp_data   = Pdfs([path_S15+'stats/pdfs/pdf76000.LnEnstrophyW_iW_i'],path_S15+'y.dat')
S20_vortpdf_interp_data   = Pdfs([path_S20_2+'stats/pdfs/pdf75000.LnEnstrophyW_iW_i'],path_S20_2+'y.dat')
S25_vortpdf_interp_data   = Pdfs([path_S25+'stats/pdfs/pdf110000.LnEnstrophyW_iW_i'],path_S25+'y.dat')


# In[7]:


# Interpolate pdfs in x-direction

S0_vortpdf_interp = np.zeros((len(vortlist_S0),S0.y_len,vortpdf_S0.nb))
for n in range(len(vortlist_S0)):
    for j in range(S0.y_len):
        S0_vortpdf_interp[n,j,:] = np.interp(S0_vortpdf_interp_data.xy[0,0,j,:],vortpdf_S0.xy[0,n,j,:],vortpdf_S0.pdf[n,j,:-2])

S05_vortpdf_interp = np.zeros((len(vortlist_S05),S05.y_len,vortpdf_S05.nb))
for n in range(len(vortlist_S05)):
    for j in range(S05.y_len):
        S05_vortpdf_interp[n,j,:] = np.interp(S05_vortpdf_interp_data.xy[0,0,j,:],vortpdf_S05.xy[0,n,j,:],vortpdf_S05.pdf[n,j,:-2])

S10_vortpdf_interp = np.zeros((len(vortlist_S10),S10.y_len,vortpdf_S10.nb))
for n in range(len(vortlist_S10)):
    for j in range(S10.y_len):
        S10_vortpdf_interp[n,j,:] = np.interp(S10_vortpdf_interp_data.xy[0,0,j,:],vortpdf_S10.xy[0,n,j,:],vortpdf_S10.pdf[n,j,:-2])

S15_vortpdf_interp = np.zeros((len(vortlist_S15),S15.y_len,vortpdf_S15.nb))
for n in range(len(vortlist_S15)):
    for j in range(S15.y_len):
        S15_vortpdf_interp[n,j,:] = np.interp(S15_vortpdf_interp_data.xy[0,0,j,:],vortpdf_S15.xy[0,n,j,:],vortpdf_S15.pdf[n,j,:-2])
        
S20_vortpdf_interp_1 = np.zeros((len(vortlist_S20_1),S20_1.y_len,vortpdf_S20_1.nb))
for n in range(len(vortlist_S20_1)):
    for j in range(S20_1.y_len):
        S20_vortpdf_interp_1[n,j,:] = np.interp(S20_vortpdf_interp_data.xy[0,0,j,:],vortpdf_S20_1.xy[0,n,j,:],vortpdf_S20_1.pdf[n,j,:-2])

S20_vortpdf_interp_2 = np.zeros((len(vortlist_S20_2),S20_2.y_len,vortpdf_S20_2.nb))
for n in range(len(vortlist_S20_2)):
    for j in range(S20_2.y_len):
        S20_vortpdf_interp_2[n,j,:] = np.interp(S20_vortpdf_interp_data.xy[0,0,j,:],vortpdf_S20_2.xy[0,n,j,:],vortpdf_S20_2.pdf[n,j,:-2])

S20_vortpdf_interp = np.concatenate((S20_vortpdf_interp_1,S20_vortpdf_interp_2),axis=0)

S25_vortpdf_interp = np.zeros((len(vortlist_S25),S25.y_len,vortpdf_S25.nb))
for n in range(len(vortlist_S25)):
    for j in range(S25.y_len):
        S25_vortpdf_interp[n,j,:] = np.interp(S25_vortpdf_interp_data.xy[0,0,j,:],vortpdf_S25.xy[0,n,j,:],vortpdf_S25.pdf[n,j,:-2])


# In[8]:


# Running mean of pdfs

S0_vortpdf_interp_runmean = np.zeros((np.ma.size(S0_vortpdf_interp,0)-2,S0.y_len,vortpdf_S0.nb))
for n in range(1,np.ma.size(S0_vortpdf_interp,0)-1):
    S0_vortpdf_interp_runmean[n-1,:,:] = np.mean(S0_vortpdf_interp[n-1:n+2,:,:],axis=0)

S05_vortpdf_interp_runmean = np.zeros((np.ma.size(S05_vortpdf_interp,0)-2,S05.y_len,vortpdf_S05.nb))
for n in range(1,np.ma.size(S05_vortpdf_interp,0)-1):
    S05_vortpdf_interp_runmean[n-1,:,:] = np.mean(S05_vortpdf_interp[n-1:n+2,:,:],axis=0)

S10_vortpdf_interp_runmean = np.zeros((np.ma.size(S10_vortpdf_interp,0)-2,S10.y_len,vortpdf_S10.nb))
for n in range(1,np.ma.size(S10_vortpdf_interp,0)-1):
    S10_vortpdf_interp_runmean[n-1,:,:] = np.mean(S10_vortpdf_interp[n-1:n+2,:,:],axis=0)

S15_vortpdf_interp_runmean = np.zeros((np.ma.size(S15_vortpdf_interp,0)-2,S15.y_len,vortpdf_S15.nb))
for n in range(1,np.ma.size(S15_vortpdf_interp,0)-1):
    S15_vortpdf_interp_runmean[n-1,:,:] = np.mean(S15_vortpdf_interp[n-1:n+2,:,:],axis=0)

S20_vortpdf_interp_runmean = np.zeros((np.ma.size(S20_vortpdf_interp,0)-2,S20_1.y_len,vortpdf_S20_1.nb))
for n in range(1,np.ma.size(S20_vortpdf_interp,0)-1):
    S20_vortpdf_interp_runmean[n-1,:,:] = np.mean(S20_vortpdf_interp[n-1:n+2,:,:],axis=0)

S25_vortpdf_interp_runmean = np.zeros((np.ma.size(S25_vortpdf_interp,0)-2,S25.y_len,vortpdf_S25.nb))
for n in range(1,np.ma.size(S25_vortpdf_interp,0)-1):
    S25_vortpdf_interp_runmean[n-1,:,:] = np.mean(S25_vortpdf_interp[n-1:n+2,:,:],axis=0)


# In[9]:


# Find mean of pdf at each height

meanvort_S0 = np.zeros((np.ma.size(S0_vortpdf_interp_runmean,0),S0.y_len))
meanprob_vort_S0 = np.zeros((np.ma.size(S0_vortpdf_interp_runmean,0),S0.y_len))
for t in range(0,np.ma.size(S0_vortpdf_interp_runmean,0)):
    for j in range(0,S0.y_len):
        meanvort_S0[t,j] = np.average(S0_vortpdf_interp_data.xy[0,0,j,:],weights=S0_vortpdf_interp_runmean[t,j,:])
        meanprob_vort_S0[t,j] = S0_vortpdf_interp_runmean[t,j,np.argmin( np.abs(S0_vortpdf_interp_data.xy[0,0,j,:] -meanvort_S0[t,j]) )]
meanvort_S0 = np.log10(np.exp(meanvort_S0)/(ceps*B0/nu))

meanvort_S05 = np.zeros((np.ma.size(S05_vortpdf_interp_runmean,0),S05.y_len))
meanprob_vort_S05 = np.zeros((np.ma.size(S05_vortpdf_interp_runmean,0),S05.y_len))
for t in range(0,np.ma.size(S05_vortpdf_interp_runmean,0)):
    for j in range(0,S05.y_len):
        meanvort_S05[t,j] = np.average(S05_vortpdf_interp_data.xy[0,0,j,:],weights=S05_vortpdf_interp_runmean[t,j,:])
        meanprob_vort_S05[t,j] = S05_vortpdf_interp_runmean[t,j,np.argmin( np.abs(S05_vortpdf_interp_data.xy[0,0,j,:] -meanvort_S05[t,j]) )]
meanvort_S05 = np.log10(np.exp(meanvort_S05)/(ceps*B0/nu))

meanvort_S10 = np.zeros((np.ma.size(S10_vortpdf_interp_runmean,0),S10.y_len))
meanprob_vort_S10 = np.zeros((np.ma.size(S10_vortpdf_interp_runmean,0),S10.y_len))
for t in range(0,np.ma.size(S10_vortpdf_interp_runmean,0)):
    for j in range(0,S10.y_len):
        meanvort_S10[t,j] = np.average(S10_vortpdf_interp_data.xy[0,0,j,:],weights=S10_vortpdf_interp_runmean[t,j,:])
        meanprob_vort_S10[t,j] = S10_vortpdf_interp_runmean[t,j,np.argmin( np.abs(S10_vortpdf_interp_data.xy[0,0,j,:] -meanvort_S10[t,j]) )]
meanvort_S10 = np.log10(np.exp(meanvort_S10)/(ceps*B0/nu))

meanvort_S15 = np.zeros((np.ma.size(S15_vortpdf_interp_runmean,0),S15.y_len))
meanprob_vort_S15 = np.zeros((np.ma.size(S15_vortpdf_interp_runmean,0),S15.y_len))
for t in range(0,np.ma.size(S15_vortpdf_interp_runmean,0)):
    for j in range(0,S15.y_len):
        meanvort_S15[t,j] = np.average(S15_vortpdf_interp_data.xy[0,0,j,:],weights=S15_vortpdf_interp_runmean[t,j,:])
        meanprob_vort_S15[t,j] = S15_vortpdf_interp_runmean[t,j,np.argmin( np.abs(S15_vortpdf_interp_data.xy[0,0,j,:] -meanvort_S15[t,j]) )]
meanvort_S15 = np.log10(np.exp(meanvort_S15)/(ceps*B0/nu))

meanvort_S20 = np.zeros((np.ma.size(S20_vortpdf_interp_runmean,0),S20_1.y_len))
meanprob_vort_S20 = np.zeros((np.ma.size(S20_vortpdf_interp_runmean,0),S20_1.y_len))
for t in range(0,np.ma.size(S20_vortpdf_interp_runmean,0)):
    for j in range(0,S20_1.y_len):
        meanvort_S20[t,j] = np.average(S20_vortpdf_interp_data.xy[0,0,j,:],weights=S20_vortpdf_interp_runmean[t,j,:])
        meanprob_vort_S20[t,j] = S20_vortpdf_interp_runmean[t,j,np.argmin( np.abs(S20_vortpdf_interp_data.xy[0,0,j,:] -meanvort_S20[t,j]) )]
meanvort_S20 = np.log10(np.exp(meanvort_S20)/(ceps*B0/nu))

meanvort_S25 = np.zeros((np.ma.size(S25_vortpdf_interp_runmean,0),S25.y_len))
meanprob_vort_S25 = np.zeros((np.ma.size(S25_vortpdf_interp_runmean,0),S25.y_len))
for t in range(0,np.ma.size(S25_vortpdf_interp_runmean,0)):
    for j in range(0,S25.y_len):
        meanvort_S25[t,j] = np.average(S25_vortpdf_interp_data.xy[0,0,j,:],weights=S25_vortpdf_interp_runmean[t,j,:])
        meanprob_vort_S25[t,j] = S25_vortpdf_interp_runmean[t,j,np.argmin( np.abs(S25_vortpdf_interp_data.xy[0,0,j,:] -meanvort_S25[t,j]) )]
meanvort_S25 = np.log10(np.exp(meanvort_S25)/(ceps*B0/nu))


# In[10]:


# Find saddle as point where meanprob has a minimum

y_vort_S0_saddle = np.zeros(np.ma.size(S0_vortpdf_interp_runmean,0))
meanvort_S0_saddle = np.zeros(np.ma.size(S0_vortpdf_interp_runmean,0))
for t in range(np.ma.size(S0_vortpdf_interp_runmean,0)):
    j_min = np.argmin(np.abs(S0.y -0.5 *z_enc_runningmean_S0[t]))
    j_max = np.argmin(np.abs(S0.y -2.0 *z_enc_runningmean_S0[t]))

    y_vort_S0_saddle[t] = S0.y[np.argmin(meanprob_vort_S0[t,j_min:j_max])+j_min]
    meanvort_S0_saddle[t] = meanvort_S0[t,np.argmin(np.abs(y_vort_S0_saddle[t]-S0.y))]
    y_vort_S0_saddle[t] = y_vort_S0_saddle[t]/z_enc_runningmean_S0[t]

y_vort_S05_saddle = np.zeros(np.ma.size(S05_vortpdf_interp_runmean,0))
meanvort_S05_saddle = np.zeros(np.ma.size(S05_vortpdf_interp_runmean,0))
for t in range(np.ma.size(S05_vortpdf_interp_runmean,0)):
    j_min = np.argmin(np.abs(S05.y -0.5 *z_enc_runningmean_S05[t]))
    j_max = np.argmin(np.abs(S05.y -2.0 *z_enc_runningmean_S05[t]))

    y_vort_S05_saddle[t] = S05.y[np.argmin(meanprob_vort_S05[t,j_min:j_max])+j_min]
    meanvort_S05_saddle[t] = meanvort_S05[t,np.argmin(np.abs(y_vort_S05_saddle[t]-S05.y))]
    y_vort_S05_saddle[t] = y_vort_S05_saddle[t]/z_enc_runningmean_S05[t]

y_vort_S10_saddle = np.zeros(np.ma.size(S10_vortpdf_interp_runmean,0))
meanvort_S10_saddle = np.zeros(np.ma.size(S10_vortpdf_interp_runmean,0))
for t in range(np.ma.size(S10_vortpdf_interp_runmean,0)):
    j_min = np.argmin(np.abs(S10.y -0.5 *z_enc_runningmean_S10[t]))
    j_max = np.argmin(np.abs(S10.y -2.0 *z_enc_runningmean_S10[t]))

    y_vort_S10_saddle[t] = S10.y[np.argmin(meanprob_vort_S10[t,j_min:j_max])+j_min]
    meanvort_S10_saddle[t] = meanvort_S10[t,np.argmin(np.abs(y_vort_S10_saddle[t]-S10.y))]
    y_vort_S10_saddle[t] = y_vort_S10_saddle[t]/z_enc_runningmean_S10[t]

y_vort_S15_saddle = np.zeros(np.ma.size(S15_vortpdf_interp_runmean,0))
meanvort_S15_saddle = np.zeros(np.ma.size(S15_vortpdf_interp_runmean,0))
for t in range(np.ma.size(S15_vortpdf_interp_runmean,0)):
    j_min = np.argmin(np.abs(S15.y -0.5 *z_enc_runningmean_S15[t]))
    j_max = np.argmin(np.abs(S15.y -2.0 *z_enc_runningmean_S15[t]))

    y_vort_S15_saddle[t] = S15.y[np.argmin(meanprob_vort_S15[t,j_min:j_max])+j_min]
    meanvort_S15_saddle[t] = meanvort_S15[t,np.argmin(np.abs(y_vort_S15_saddle[t]-S15.y))]
    y_vort_S15_saddle[t] = y_vort_S15_saddle[t]/z_enc_runningmean_S15[t]

y_vort_S20_saddle = np.zeros(np.ma.size(S20_vortpdf_interp_runmean,0))
meanvort_S20_saddle = np.zeros(np.ma.size(S20_vortpdf_interp_runmean,0))
for t in range(np.ma.size(S20_vortpdf_interp_runmean,0)):
    j_min = np.argmin(np.abs(S20_1.y -0.5 *z_enc_runningmean_S20[t]))
    j_max = np.argmin(np.abs(S20_1.y -2.0 *z_enc_runningmean_S20[t]))

    y_vort_S20_saddle[t] = S20_1.y[np.argmin(meanprob_vort_S20[t,j_min:j_max])+j_min]
    meanvort_S20_saddle[t] = meanvort_S20[t,np.argmin(np.abs(y_vort_S20_saddle[t]-S20_1.y))]
    y_vort_S20_saddle[t] = y_vort_S20_saddle[t]/z_enc_runningmean_S20[t]

y_vort_S25_saddle = np.zeros(np.ma.size(S25_vortpdf_interp_runmean,0))
meanvort_S25_saddle = np.zeros(np.ma.size(S25_vortpdf_interp_runmean,0))
for t in range(np.ma.size(S25_vortpdf_interp_runmean,0)):
    j_min = np.argmin(np.abs(S25.y -0.5 *z_enc_runningmean_S25[t]))
    j_max = np.argmin(np.abs(S25.y -2.0 *z_enc_runningmean_S25[t]))

    y_vort_S25_saddle[t] = S25.y[np.argmin(meanprob_vort_S25[t,j_min:j_max])+j_min]
    meanvort_S25_saddle[t] = meanvort_S25[t,np.argmin(np.abs(y_vort_S25_saddle[t]-S25.y))]
    y_vort_S25_saddle[t] = y_vort_S25_saddle[t]/z_enc_runningmean_S25[t]


# In[11]:


# concatenate over time

time = np.concatenate((S20_1.z_enc/L0,S20_2.z_enc/L0))
z_if = np.concatenate((S20_1.z_if/S20_1.z_enc,S20_2.z_if/S20_2.z_enc))
z_ig = np.concatenate((S20_1.z_ig/S20_1.z_enc,S20_2.z_ig/S20_2.z_enc))


# In[13]:


# Plot

#blues = matplotlib.cm.get_cmap('Blues')

colors = []
for value in [0.1,0.25,0.40,0.55,0.7,0.85]:
    colors.append(matplotlib.cm.get_cmap('plasma_r')(value))

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax1.set_ylim(1,1.6)
ax1.set_xlim(15,30)
ax2.set_ylim(-2,0)
ax2.set_xlim(15,30)
ax1.plot(S0.z_enc[1:-1]/L0, y_vort_S0_saddle,   c=colors[0])
ax1.plot(S05.z_enc[1:-1]/L0,y_vort_S05_saddle,  c=colors[1])
ax1.plot(S10.z_enc[1:-1]/L0,y_vort_S10_saddle,  c=colors[2])
ax1.plot(S15.z_enc[1:-1]/L0,y_vort_S15_saddle,  c=colors[3])
ax1.plot(time[1:-1],        y_vort_S20_saddle,  c=colors[4])
ax1.plot(S25.z_enc[1:-1]/L0,y_vort_S25_saddle,  c=colors[5])
ax2.plot(S0.z_enc[1:-1]/L0, meanvort_S0_saddle,  c=colors[0], label=r'$Fr_0=0$')
ax2.plot(S05.z_enc[1:-1]/L0,meanvort_S05_saddle, c=colors[1], label=r'$Fr_0=5$')
ax2.plot(S10.z_enc[1:-1]/L0,meanvort_S10_saddle, c=colors[2], label=r'$Fr_0=10$')
ax2.plot(S15.z_enc[1:-1]/L0,meanvort_S15_saddle, c=colors[3], label=r'$Fr_0=15$')
ax2.plot(time[1:-1],        meanvort_S20_saddle, c=colors[4], label=r'$Fr_0=20$')
ax2.plot(S25.z_enc[1:-1]/L0,meanvort_S25_saddle, c=colors[5], label=r'$Fr_0=25$')
ax1.set_xlabel(r'$z_{enc}/L_0$')
ax2.set_xlabel(r'$z_{enc}/L_0$')
ax1.set_ylabel(r'$z_\omega/z_{enc}$')
ax2.set_ylabel(r'$\mathrm{log}_{10}(\omega^2/\omega_0^2)$')
ax1.set_title('(a)',fontsize=24,loc='left')
ax2.set_title('(b)',fontsize=24,loc='left')
ax2.legend(loc='best',fontsize=18,handlelength=1,borderaxespad=0.2,ncol=2,columnspacing=1)
plt.tight_layout()
plt.savefig(opath+'pdfs_saddle_point_mean.pdf',bbox_inches='tight')
plt.show()


# In[14]:


np.mean(meanvort_S0_saddle)


# In[15]:


np.mean(meanvort_S05_saddle)


# In[19]:


np.mean(meanvort_S10_saddle)


# In[20]:


np.mean(meanvort_S15_saddle)


# In[21]:


np.mean(meanvort_S20_saddle)


# In[22]:


np.mean(meanvort_S25_saddle)


# In[5]:


np.sqrt((10**(-0.46751851045102605))*(ceps*B0/nu))


# In[16]:


np.mean(y_vort_S0_saddle)


# In[ ]:




