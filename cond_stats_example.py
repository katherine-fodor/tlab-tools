#!/usr/bin/python3
# Example script for plotting conditional statistics using ReadStats.
# Data is from a convective boundary layer (CBL) driven by a surface buoyancy flux, B0,
# growing into an overlying stratification with buoyancy gradient, N^2.
# Example is for plotting the intermittency factor and the buoyancy flux within turbulent
# regions for a shear-free CBL and a sheared CBL.

#####################################################################
# Import libraries, set parameters for figures.
# Set an output path for plots.

from ReadStats import Statistics, Conditional_Stats
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.lines as mlines
from matplotlib import rc

rc('text', usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font', family='serif')
rc('font', size=24)
rc('axes', linewidth=1.5)
rc('lines', linewidth=2)

opath = '/scratch/local1/m300551/ForKatherine/plots/3D/Re042/'

#######################################################################
# Set some useful constants.

nu = 1./25000.
B0 = 0.005
N = np.sqrt(3)
L0 = (B0/N**3)**0.5

#######################################################################
# Calculate a running mean of a time series with a specified window size.
# Window size is number of entries on either side of the entry being averaged.
# Leaves out entries at the beginning and end of the time series such that the
# window size is always the same, but the resulting time series is shorter than
# the original.

def runningmean(timeseries,window):
    nt = len(timeseries)
    outseries = np.zeros(nt-(window*2))
    for n in range(window,nt-window):
        outseries[n-window] = np.mean(timeseries[n-window:n+window+1])
    return outseries

#######################################################################
# Set the paths to the required data

path_1 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/2560x576x2560/'
path_2 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/2560x704x2560/'
path_3 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/2560x896x2560/'

path_S20 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re042/3072x960x4608-S20/'

path_vort = 'stats/gate-vorticity/gate-2-08/'
path_vort_S20 = 'stats/gate-vorticity/gate-1-24/'

# Use the class Statistics to read in conventional statistics from Tlab.

NS42_1 = Statistics(path_1+'stats/pdftimes/avg20500-53000.nc')
NS42_2 = Statistics(path_2+'stats/pdftimes/avg60000-74500.nc')
NS42_3 = Statistics(path_3+'stats/pdftimes/avg83000-127500.nc')

S20 = Statistics(path_S20+'stats/pdftimes/avg42000-148000.nc')

# Use the class Conditional_Stats to read in conditional statistics from Tlab.

NS42_vort_int_1 = Conditional_Stats(path_1+path_vort+'int20500-53000.nc',path_1+path_vort+'Partition1/cavg20500-53000.nc',path_1+path_vort+'Partition2/cavg20500-53000.nc')
NS42_vort_int_2 = Conditional_Stats(path_2+path_vort+'int60000-74500.nc',path_2+path_vort+'Partition1/cavg60000-74500.nc',path_2+path_vort+'Partition2/cavg60000-74500.nc')
NS42_vort_int_3 = Conditional_Stats(path_3+path_vort+'int83000-127500.nc',path_3+path_vort+'Partition1/cavg83000-127500.nc',path_3+path_vort+'Partition2/cavg83000-127500.nc')

S20_vort_int = Conditional_Stats(path_S20+path_vort_S20+'int42000-148000.nc',path_S20+path_vort_S20+'Partition1/cavg42000-148000.nc',path_S20+path_vort_S20+'Partition2/cavg42000-148000.nc')

# Call the intermittency factor.

NS42_vort_turbareafrac_zif_1 = [NS42_vort_int_1.int2[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_turbareafrac_zif_2 = [NS42_vort_int_2.int2[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_turbareafrac_zif_3 = [NS42_vort_int_3.int2[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_turbareafrac_zif = NS42_vort_turbareafrac_zif_1 + NS42_vort_turbareafrac_zif_2 +  NS42_vort_turbareafrac_zif_3

S20_vort_turbareafrac_zif = [S20_vort_int.int2[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# Call some variables in Partition 1 (in this case, statistics conditioned to non-turbulent regions).

NS42_vort_p1_s1_mean_zif_1 = [NS42_vort_int_1.P1S1Mom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_s1_mean_zif_2 = [NS42_vort_int_2.P1S1Mom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_s1_mean_zif_3 = [NS42_vort_int_3.P1S1Mom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_s1_mean_zif = NS42_vort_p1_s1_mean_zif_1 + NS42_vort_p1_s1_mean_zif_2 + NS42_vort_p1_s1_mean_zif_3

NS42_vort_p1_w_mean_zif_1 = [NS42_vort_int_1.P1VMom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_w_mean_zif_2 = [NS42_vort_int_2.P1VMom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_w_mean_zif_3 = [NS42_vort_int_3.P1VMom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_w_mean_zif = NS42_vort_p1_w_mean_zif_1 + NS42_vort_p1_w_mean_zif_2 + NS42_vort_p1_w_mean_zif_3

NS42_vort_p1_v1_zif_1 = [NS42_vort_int_1.P1v1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p1_v1_zif_2 = [NS42_vort_int_2.P1v1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p1_v1_zif_3 = [NS42_vort_int_3.P1v1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p1_v1_zif = NS42_vort_p1_v1_zif_1 + NS42_vort_p1_v1_zif_2 + NS42_vort_p1_v1_zif_3

S20_vort_p1_s1_mean_zif = [S20_vort_int.P1S1Mom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_w_mean_zif = [S20_vort_int.P1VMom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p1_v1_zif = [S20_vort_int.P1v1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# Call some variables in Partition 2 (in this case, statistics conditioned to turbulent regions).

NS42_vort_p2_s1_mean_zif_1 = [NS42_vort_int_1.P2S1Mom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_s1_mean_zif_2 = [NS42_vort_int_2.P2S1Mom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_s1_mean_zif_3 = [NS42_vort_int_3.P2S1Mom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_s1_mean_zif = NS42_vort_p2_s1_mean_zif_1 + NS42_vort_p2_s1_mean_zif_2 + NS42_vort_p2_s1_mean_zif_3

NS42_vort_p2_w_mean_zif_1 = [NS42_vort_int_1.P2VMom1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_w_mean_zif_2 = [NS42_vort_int_2.P2VMom1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_w_mean_zif_3 = [NS42_vort_int_3.P2VMom1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_w_mean_zif = NS42_vort_p2_w_mean_zif_1 + NS42_vort_p2_w_mean_zif_2 + NS42_vort_p2_w_mean_zif_3

NS42_vort_p2_v1_zif_1 = [NS42_vort_int_1.P2v1[n,NS42_1.z_if_arg[n]] for n in range(0,NS42_1.t_len)]
NS42_vort_p2_v1_zif_2 = [NS42_vort_int_2.P2v1[n,NS42_2.z_if_arg[n]] for n in range(0,NS42_2.t_len)]
NS42_vort_p2_v1_zif_3 = [NS42_vort_int_3.P2v1[n,NS42_3.z_if_arg[n]] for n in range(0,NS42_3.t_len)]

NS42_vort_p2_v1_zif = NS42_vort_p2_v1_zif_1 + NS42_vort_p2_v1_zif_2 + NS42_vort_p2_v1_zif_3

S20_vort_p2_s1_mean_zif = [S20_vort_int.P2S1Mom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_w_mean_zif = [S20_vort_int.P2VMom1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]
S20_vort_p2_v1_zif = [S20_vort_int.P2v1[n,S20.z_if_arg[n]] for n in range(0,S20.t_len)]

# Call some useful variables from Statistics.

time_1 = [NS42_1.z_enc[n]/L0 for n in range(0,NS42_1.t_len)]
time_2 = [NS42_2.z_enc[n]/L0 for n in range(0,NS42_2.t_len)]
time_3 = [NS42_3.z_enc[n]/L0 for n in range(0,NS42_3.t_len)]

time = time_1 + time_2 + time_3

z_if = np.concatenate([NS42_1.z_if,NS42_2.z_if,NS42_3.z_if])
z_ig = np.concatenate([NS42_1.z_ig,NS42_2.z_ig,NS42_3.z_ig])
z_enc = np.concatenate([NS42_1.z_enc,NS42_2.z_enc,NS42_3.z_enc])
z_is = np.concatenate([NS42_1.z_is,NS42_2.z_is,NS42_3.z_is])
L_Oz_zif = np.concatenate([NS42_1.L_Oz_zif,NS42_2.L_Oz_zif,NS42_3.L_Oz_zif])
w_enc = (B0*z_enc)**(1./3.)

#######################################################################
# Plot

blues = matplotlib.cm.get_cmap('Blues')

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax1.set_ylim(0,1.6)
ax2.set_ylim(0.4,1)
ax2.set_xlim(15,30)
ax1.plot(np.mean(NS42_vort_int_1.int2[-4:-1,:],axis=0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),label=r'$Fr_0=0$')
ax1.plot(np.mean(S20_vort_int.int2[4:7,:],axis=0),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),label=r'$Fr_0=20$')
ax2.plot(time[1:-1],runningmean(NS42_vort_turbareafrac_zif,1),c=blues(0.5))
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(S20_vort_turbareafrac_zif,1),c=blues(0.9))
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$a_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$(a_\mathrm{T})_{z_{i,f}}$')
ax2.set_title(r'(b)',fontsize=20,loc='left')
ax1.legend(loc='lower left',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'area_frac_height_time_S20_S0_zif_vort.pdf',bbox_inches='tight')
plt.show()

f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True)
ax2.grid(True)
ax1.set_ylim(1,1.4)
ax1.set_xlim(-0.4,0.2)
ax2.set_ylim(0,0.4)
ax2.set_xlim(15,30)
ax1.plot(np.mean(NS42_vort_int_1.P2v1[-4:-1,:]-NS42_vort_int_1.P2S1Mom1[-4:-1,:]*NS42_vort_int_1.P2VMom1[-4:-1,:],axis=0)/B0,NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5))
ax1.plot(np.mean(S20_vort_int.P2v1[4:7,:]-S20_vort_int.P2S1Mom1[4:7,:]*S20_vort_int.P2VMom1[4:7,:],axis=0)/B0,S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9))
ax2.plot(time[1:-1],-runningmean(np.array(NS42_vort_p2_v1_zif)-np.array(NS42_vort_p2_s1_mean_zif)*np.array(NS42_vort_p2_w_mean_zif),1)/B0,c=blues(0.5),label=r'$Fr_0=0$')
ax2.plot(S20.z_enc[1:-1]/L0,-runningmean(np.array(S20_vort_p2_v1_zif)-np.array(S20_vort_p2_s1_mean_zif)*np.array(S20_vort_p2_w_mean_zif),1)/B0,c=blues(0.9),label=r'$Fr_0=20$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$\langle b^\prime w^\prime\rangle_\mathrm{T}/B_0$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$-(\langle b^\prime w^\prime \rangle_\mathrm{T})_{z_{i,f}}/B_0$')
ax2.set_title(r'(b)',fontsize=20,loc='left')
ax2.legend(loc='best',fontsize=20)
plt.tight_layout()
plt.savefig(opath+'s1_vflux_turb_height_time_S20_S0_vort.pdf')
plt.show()



