#!/usr/bin/python3
# Example script for plotting PDFs using ReadStats.
# Data is from a convective boundary layer (CBL) driven by a surface buoyancy flux, B0,
# growing into an overlying stratification with buoyancy gradient, N^2.
# Example is for plotting the enstrophy PDF for a shear-free CBL and a sheared CBL.

##########################################################################
# Import libraries, set parameters for figures.
# Set an output path for plots and an import path for data.

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

opath = '/scratch/local1/m300551/ForKatherine/plots/3D/Re042/'
path = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/'
##########################################################################
# Set some useful constants.

nu = 1./25000.
B0 = 0.005
N = np.sqrt(3)
L0 = (B0/N**3)**0.5
ceps = 0.1

##########################################################################
# Use the class Statistics to read in conventional statistics from Tlab.

NS42 = Statistics(path+'2560x576x2560/stats/pdftimes/avg36500-47500.nc')
NS42_3 = Statistics(path+'2560x896x2560/stats/pdftimes/avg83000-127500.nc')
S20 = Statistics(path+'3072x960x4608-S20/stats/pdftimes/avg66000-84000.nc')

############################################################################
# Use the class Pdfs to read in the PDF data from Tlab.

NS42_vortpdf  = Pdfs([path+'2560x576x2560/stats/pdfs/pdf36500.LnEnstrophyW_iW_i',path+'2560x576x2560/stats/pdfs/pdf41500.LnEnstrophyW_iW_i',path+'2560x576x2560/stats/pdfs/pdf47500.LnEnstrophyW_iW_i'],path+'2560x576x2560/y.dat')
S20_vortpdf = Pdfs([path+'3072x960x4608-S20/stats/pdfs/pdf66000.LnEnstrophyW_iW_i',path+'3072x960x4608-S20/stats/pdfs/pdf75000.LnEnstrophyW_iW_i',path+'3072x960x4608-S20/stats/pdfs/pdf84000.LnEnstrophyW_iW_i'],path+'3072x960x4608-S20/y.dat')

# Create grid on which to interpolate PDFs

NS42_vortpdf_interp_data = Pdfs([path+'2560x896x2560/stats/pdfs/pdf102500.LnEnstrophyW_iW_i'],path+'2560x896x2560/y.dat')
S20_vortpdf_interp_data = Pdfs([path+'3072x960x4608-S20/stats/pdfs/pdf75000.LnEnstrophyW_iW_i'],path+'3072x960x4608-S20/y.dat')

# Interpolate PDFs in y-direction (only needed if using data with different vertical grids)

NS42_vortpdf_interp_y = np.zeros((3,NS42_3.y_len,NS42_vortpdf.nb+2))
for n in range(3):
    for i in range(NS42_vortpdf.nb+2):
        NS42_vortpdf_interp_y[n,:,i] = np.interp(NS42_3.y,NS42.y,NS42_vortpdf.pdf[n,:-1,i])

# Interpolate PDFs in x-direction

NS42_vortpdf_interp = np.zeros((3,NS42_3.y_len,NS42_vortpdf.nb))
for n in range(3):
    for j in range(NS42_3.y_len):
        NS42_vortpdf_interp[n,j,:] = np.interp(NS42_vortpdf_interp_data.xy[0,0,j,:],np.linspace(NS42_vortpdf_interp_y[n,j,NS42_vortpdf.nb],NS42_vortpdf_interp_y[n,j,NS42_vortpdf.nb+1],num=NS42_vortpdf.nb),NS42_vortpdf_interp_y[n,j,:-2])


S20_vortpdf_interp = np.zeros((3,S20.y_len,S20_vortpdf.nb))
for n in range(3):
    for j in range(S20.y_len):
        S20_vortpdf_interp[n,j,:] = np.interp(S20_vortpdf_interp_data.xy[0,0,j,:],S20_vortpdf.xy[0,n,j,:],S20_vortpdf.pdf[n,j,:-2])

# Calculate the mean of the PDFs

NS42_vortpdf_interp_runmean = np.mean(NS42_vortpdf_interp,axis=0)
S20_vortpdf_interp_runmean = np.mean(S20_vortpdf_interp,axis=0)

# Find where PDF has a maximum at each height

maxvort_NS42 = np.zeros(NS42_3.y_len)
maxprob_vort_NS42 = np.zeros(NS42_3.y_len)
for i in range(0,NS42_3.y_len):
    maxvort_NS42[i] = NS42_vortpdf_interp_data.xy[0,0,i,np.argmax(NS42_vortpdf_interp_runmean[i,:])]
    maxprob_vort_NS42[i] = np.max(NS42_vortpdf_interp_runmean[i,:])
maxvort_NS42 = np.log10(np.exp(maxvort_NS42)/(ceps*B0/nu))

maxvort_S20 = np.zeros(S20.y_len)
maxprob_vort_S20 = np.zeros(S20.y_len)
for i in range(0,S20.y_len):
    maxvort_S20[i] = S20_vortpdf_interp_data.xy[0,0,i,np.argmax(S20_vortpdf_interp_runmean[i,:])]
    maxprob_vort_S20[i] = np.max(S20_vortpdf_interp_runmean[i,:])
maxvort_S20 = np.log10(np.exp(maxvort_S20)/(ceps*B0/nu))

# Find jump in maxvort

for i in range(NS42_3.y_len):
    if np.abs(maxvort_NS42[i+1])-np.abs(maxvort_NS42[i]) > 0.2:
        maxit_vort_NS42 = i+1
        break
  
for i in range(S20.y_len):
    if np.abs(maxvort_S20[i+1])-np.abs(maxvort_S20[i]) > 0.5:
        maxit_vort_S20 = i+1
        break
    
# Find saddle point as point where maxprob has a minimum

y_vort_NS42_saddle = NS42_3.y[np.argmin(maxprob_vort_NS42[:-100])]
maxvort_NS42_saddle = maxvort_NS42[np.argmin(np.abs(y_vort_NS42_saddle-NS42_3.y))]
y_vort_NS42_saddle = y_vort_NS42_saddle/np.mean(NS42.z_enc)

y_vort_S20_saddle = S20.y[np.argmin(maxprob_vort_S20)]
maxvort_S20_saddle = maxvort_S20[np.argmin(np.abs(y_vort_S20_saddle-S20.y))]
y_vort_S20_saddle = y_vort_S20_saddle/np.mean(S20.z_enc)

# Normalise y-axis

NS42_vortpdf_y_mean = NS42_vortpdf_interp_data.xy[1,0,:,:]/np.mean(NS42.z_enc)
S20_42_vortpdf_y_mean = S20_vortpdf_interp_data.xy[1,0,:,:]/np.mean(S20.z_enc)

# Normalise x-axis

NS42_vortpdf_x_mean = np.log10(np.exp(NS42_vortpdf_interp_data.xy[0,0,:,:])/(ceps*B0/nu_))
S20_42_vortpdf_x_mean = np.log10(np.exp(S20_vortpdf_interp_data.xy[0,0,:,:])/(ceps*B0/nu))

#####################################################################
# Plot

f, (ax1,ax2) = plt.subplots(1,2,sharex='all',sharey='all',figsize=(10,5))
ax1.tick_params(bottom=True,top=True,left=True,right=True)
ax2.tick_params(bottom=True,top=True,left=True,right=True)
ax1.set_xlim(-4,2) 
ax1.set_ylim(0,1.6)
ax1.set_xticks([-4,-3,-2,-1,0,1,2])
ax1.set_yticks([0,0.25,0.5,0.75,1,1.25,1.5])
cs1 = ax1.contourf(NS42_vortpdf_x_mean,NS42_vortpdf_y_mean,NS42_vortpdf_interp_runmean,cmap='viridis',levels=np.linspace(0,0.4,9))
ax1.plot(maxvort_NS42[:maxit_vort_NS42],NS42_3.y[:maxit_vort_NS42]/np.mean(NS42.z_enc),'k',lw=1)
ax1.plot(maxvort_NS42[maxit_vort_NS42:],NS42_3.y[maxit_vort_NS42:]/np.mean(NS42.z_enc),'k',lw=1)
ax1.scatter((maxvort_NS42[maxit_vort_NS42-1]+maxvort_NS42[maxit_vort_NS42])/2,NS42_3.y[maxit_vort_NS42]/np.mean(NS42.z_enc),100,color='k',marker='*')
ax1.axhline(np.mean(NS42.z_ig/NS42.z_enc),0,0.05,color='C1',linewidth=2)
ax1.axhline(np.mean(NS42.z_is/NS42.z_enc),0,0.05,color='C1',linewidth=2)
ax1.axhline(np.mean(NS42.z_if/NS42.z_enc),0,0.05,color='C1',linewidth=2)
cs2 = ax2.contourf(S20_42_vortpdf_x_mean,S20_42_vortpdf_y_mean,S20_vortpdf_interp_runmean,cmap='viridis',levels=np.linspace(0,0.4,9))
ax2.plot(maxvort_S20[:maxit_vort_S20],S20_42.y[:maxit_vort_S20]/np.mean(S20_42.z_enc),'k',lw=1)
ax2.plot(maxvort_S20[maxit_vort_S20:],S20_42.y[maxit_vort_S20:]/np.mean(S20_42.z_enc),'k',lw=1)
ax2.scatter((maxvort_S20[maxit_vort_S20-1]+maxvort_S20[maxit_vort_S20])/2,S20_42.y[maxit_vort_S20]/np.mean(S20_42.z_enc),100,color='k',marker='*')
ax2.axhline(np.mean(S20_42.z_ig/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(S20_42.z_is/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax2.axhline(np.mean(S20_42.z_if/S20_42.z_enc),0,0.05,color='C1',linewidth=2)
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax2.set_xlabel(r'$\log_{10}(\omega^2/\omega_0^2)$')
ax1.set_title(r'(a)$Fr_0=0$',fontsize=20,loc='left')
ax2.set_title(r'(b)$Fr_0=20$',fontsize=20,loc='left')
cbar_ax = f.add_axes([0.3,0.1,0.5,0.03])
cbar = f.colorbar(cs1,cax=cbar_ax,orientation='horizontal')
plt.tight_layout(rect=[0,0.1,1,1],h_pad=2)
plt.savefig(opath_42+'pdfs_vort_S20_S0_timeavg_interpxy.pdf',bbox_inches='tight')
plt.show()


