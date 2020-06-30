#####################################################################
# Modules

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
rc('axes', labelsize=24)
rc('lines', linewidth=2)

opath = '/scratch/local1/m300551/ForKatherine/plots/3D/Re025/Rapids/'

#######################################################################
# Constants

nu = 1./15000.
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
# Stats

path_S0 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re025/2560x512x2560/'
path_S05 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re025/1280x512x1280-S05/'
path_S10 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re025/1280x512x1280-S10/'
path_S15 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re025/1536x576x1536-S15/'
path_S20_1 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re025/S20-1536x576x1536/'
path_S20_2 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re025/S20-1536x576x2304/'
path_S25 = '/scratch/local1/m300551/ForKatherine/qCBL_3D/Re025/2560x896x2560-S25/'

path_vort_S0 = 'stats/gate-vorticity/gate-1-76/'
path_vort_S05 = 'stats/gate-vorticity/gate-2-077/'
#path_vort_S05 = 'stats/gate-vorticity/gate-1-76/'
path_vort_S10 = 'stats/gate-vorticity/gate-2-088/'
path_vort_S15 = 'stats/gate-vorticity/gate-1-056/'
path_vort_S20 = 'stats/gate-vorticity/gate-1-1/'
path_vort_S25 = 'stats/gate-vorticity/gate-1-26/'

# Conventional

S0 = Statistics(path_S0+'stats/pdftimes/avg1000-54000.nc')
S05 = Statistics(path_S05+'stats/pdftimes/avg1000-67000.nc')
S10 = Statistics(path_S10+'stats/pdftimes/avg13000-84000.nc')
S15 = Statistics(path_S15+'stats/pdftimes/avg15000-92000.nc')
S20_1 = Statistics(path_S20_1+'stats/pdftimes/avg17000-35000.nc')
S20_2 = Statistics(path_S20_2+'stats/pdftimes/avg39000-91000.nc')
S25 = Statistics(path_S25+'stats/pdftimes/avg28000-128000.nc')

S0_s1_var_zif = [S0.r2S[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]
S0_w_var_zif = [S0.Ryy[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]
S0_s1_flux_zif = [S0.Rsv[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]

S05_s1_var_zif = [S05.r2S[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]
S05_w_var_zif = [S05.Ryy[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]
S05_s1_flux_zif = [S05.Rsv[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]

S10_s1_var_zif = [S10.r2S[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]
S10_w_var_zif = [S10.Ryy[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]
S10_s1_flux_zif = [S10.Rsv[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]

S15_s1_var_zif = [S15.r2S[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]
S15_w_var_zif = [S15.Ryy[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]
S15_s1_flux_zif = [S15.Rsv[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]

S20_s1_var_zif_1 = [S20_1.r2S[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_s1_var_zif_2 = [S20_2.r2S[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_s1_var_zif = S20_s1_var_zif_1 + S20_s1_var_zif_2

S20_w_var_zif_1 = [S20_1.Ryy[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_w_var_zif_2 = [S20_2.Ryy[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_w_var_zif = S20_w_var_zif_1 + S20_w_var_zif_2

S20_s1_flux_zif_1 = [S20_1.Rsv[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_s1_flux_zif_2 = [S20_2.Rsv[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_s1_flux_zif = S20_s1_flux_zif_1 + S20_s1_flux_zif_2

S25_s1_var_zif = [S25.r2S[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]
S25_w_var_zif = [S25.Ryy[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]
S25_s1_flux_zif = [S25.Rsv[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]

## Conditional ##

S0_vort_int = Conditional_Stats(path_S0+path_vort_S0+'int1000-54000.nc',path_S0+path_vort_S0+'Partition1/cavg1000-54000.nc',path_S0+path_vort_S0+'Partition2/cavg1000-54000.nc')
S05_vort_int = Conditional_Stats(path_S05+path_vort_S05+'int1000-67000.nc',path_S05+path_vort_S05+'Partition1/cavg1000-67000.nc',path_S05+path_vort_S05+'Partition2/cavg1000-67000.nc')
S10_vort_int = Conditional_Stats(path_S10+path_vort_S10+'int13000-84000.nc',path_S10+path_vort_S10+'Partition1/cavg13000-84000.nc',path_S10+path_vort_S10+'Partition2/cavg13000-84000.nc')
S15_vort_int = Conditional_Stats(path_S15+path_vort_S15+'int15000-92000.nc',path_S15+path_vort_S15+'Partition1/cavg15000-92000.nc',path_S15+path_vort_S15+'Partition2/cavg15000-92000.nc')
S20_vort_int_1 = Conditional_Stats(path_S20_1+path_vort_S20+'int17000-35000.nc',path_S20_1+path_vort_S20+'Partition1/cavg17000-35000.nc',path_S20_1+path_vort_S20+'Partition2/cavg17000-35000.nc')
S20_vort_int_2 = Conditional_Stats(path_S20_2+path_vort_S20+'int39000-91000.nc',path_S20_2+path_vort_S20+'Partition1/cavg39000-91000.nc',path_S20_2+path_vort_S20+'Partition2/cavg39000-91000.nc')
S25_vort_int = Conditional_Stats(path_S25+path_vort_S25+'int28000-128000.nc',path_S25+path_vort_S25+'Partition1/cavg28000-128000.nc',path_S25+path_vort_S25+'Partition2/cavg28000-128000.nc')


# Vorticity #

S0_vort_turbareafrac_zif = [S0_vort_int.int2[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]
S05_vort_turbareafrac_zif = [S05_vort_int.int2[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]
S10_vort_turbareafrac_zif = [S10_vort_int.int2[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]
S15_vort_turbareafrac_zif = [S15_vort_int.int2[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]
S20_vort_turbareafrac_zif_1 = [S20_vort_int_1.int2[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_vort_turbareafrac_zif_2 = [S20_vort_int_2.int2[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_vort_turbareafrac_zif = S20_vort_turbareafrac_zif_1 + S20_vort_turbareafrac_zif_2
S25_vort_turbareafrac_zif = [S25_vort_int.int2[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]

# Non-turbulent

S0_vort_p1_s1_mean_zif = [S0_vort_int.P1S1Mom1[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]
S0_vort_p1_w_mean_zif = [S0_vort_int.P1VMom1[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]
S0_vort_p1_s1_var_zif = [S0_vort_int.P1S1Mom2[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]
S0_vort_p1_w_var_zif = [S0_vort_int.P1VMom2[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]
S0_vort_p1_v1_zif = [S0_vort_int.P1v1[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]

S05_vort_p1_s1_mean_zif = [S05_vort_int.P1S1Mom1[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]
S05_vort_p1_w_mean_zif = [S05_vort_int.P1VMom1[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]
S05_vort_p1_s1_var_zif = [S05_vort_int.P1S1Mom2[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]
S05_vort_p1_w_var_zif = [S05_vort_int.P1VMom2[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]
S05_vort_p1_v1_zif = [S05_vort_int.P1v1[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]

S10_vort_p1_s1_mean_zif = [S10_vort_int.P1S1Mom1[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]
S10_vort_p1_w_mean_zif = [S10_vort_int.P1VMom1[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]
S10_vort_p1_s1_var_zif = [S10_vort_int.P1S1Mom2[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]
S10_vort_p1_w_var_zif = [S10_vort_int.P1VMom2[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]
S10_vort_p1_v1_zif = [S10_vort_int.P1v1[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]

S15_vort_p1_s1_mean_zif = [S15_vort_int.P1S1Mom1[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]
S15_vort_p1_w_mean_zif = [S15_vort_int.P1VMom1[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]
S15_vort_p1_s1_var_zif = [S15_vort_int.P1S1Mom2[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]
S15_vort_p1_w_var_zif = [S15_vort_int.P1VMom2[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]
S15_vort_p1_v1_zif = [S15_vort_int.P1v1[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]

S20_vort_p1_s1_mean_zif_1 = [S20_vort_int_1.P1S1Mom1[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_vort_p1_s1_mean_zif_2 = [S20_vort_int_2.P1S1Mom1[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_vort_p1_s1_mean_zif = S20_vort_p1_s1_mean_zif_1 + S20_vort_p1_s1_mean_zif_2

S20_vort_p1_w_mean_zif_1 = [S20_vort_int_1.P1VMom1[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_vort_p1_w_mean_zif_2 = [S20_vort_int_2.P1VMom1[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_vort_p1_w_mean_zif = S20_vort_p1_w_mean_zif_1 + S20_vort_p1_w_mean_zif_2

S20_vort_p1_s1_var_zif_1 = [S20_vort_int_1.P1S1Mom2[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_vort_p1_s1_var_zif_2 = [S20_vort_int_2.P1S1Mom2[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_vort_p1_s1_var_zif = S20_vort_p1_s1_var_zif_1 + S20_vort_p1_s1_var_zif_2

S20_vort_p1_w_var_zif_1 = [S20_vort_int_1.P1VMom2[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_vort_p1_w_var_zif_2 = [S20_vort_int_2.P1VMom2[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_vort_p1_w_var_zif = S20_vort_p1_w_var_zif_1 + S20_vort_p1_w_var_zif_2

S20_vort_p1_v1_zif_1 = [S20_vort_int_1.P1v1[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_vort_p1_v1_zif_2 = [S20_vort_int_2.P1v1[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_vort_p1_v1_zif = S20_vort_p1_v1_zif_1 + S20_vort_p1_v1_zif_2

S25_vort_p1_s1_mean_zif = [S25_vort_int.P1S1Mom1[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]
S25_vort_p1_w_mean_zif = [S25_vort_int.P1VMom1[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]
S25_vort_p1_s1_var_zif = [S25_vort_int.P1S1Mom2[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]
S25_vort_p1_w_var_zif = [S25_vort_int.P1VMom2[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]
S25_vort_p1_v1_zif = [S25_vort_int.P1v1[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]

# Turbulent

S0_vort_p2_s1_mean_zif = [S0_vort_int.P2S1Mom1[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]
S0_vort_p2_w_mean_zif = [S0_vort_int.P2VMom1[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]
S0_vort_p2_s1_var_zif = [S0_vort_int.P2S1Mom2[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]
S0_vort_p2_w_var_zif = [S0_vort_int.P2VMom2[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]
S0_vort_p2_v1_zif = [S0_vort_int.P2v1[n,S0.z_if_arg[n]] for n in range(0,S0.t_len)]

S05_vort_p2_s1_mean_zif = [S05_vort_int.P2S1Mom1[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]
S05_vort_p2_w_mean_zif = [S05_vort_int.P2VMom1[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]
S05_vort_p2_s1_var_zif = [S05_vort_int.P2S1Mom2[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]
S05_vort_p2_w_var_zif = [S05_vort_int.P2VMom2[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]
S05_vort_p2_v1_zif = [S05_vort_int.P2v1[n,S05.z_if_arg[n]] for n in range(0,S05.t_len)]

S10_vort_p2_s1_mean_zif = [S10_vort_int.P2S1Mom1[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]
S10_vort_p2_w_mean_zif = [S10_vort_int.P2VMom1[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]
S10_vort_p2_s1_var_zif = [S10_vort_int.P2S1Mom2[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]
S10_vort_p2_w_var_zif = [S10_vort_int.P2VMom2[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]
S10_vort_p2_v1_zif = [S10_vort_int.P2v1[n,S10.z_if_arg[n]] for n in range(0,S10.t_len)]

S15_vort_p2_s1_mean_zif = [S15_vort_int.P2S1Mom1[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]
S15_vort_p2_w_mean_zif = [S15_vort_int.P2VMom1[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]
S15_vort_p2_s1_var_zif = [S15_vort_int.P2S1Mom2[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]
S15_vort_p2_w_var_zif = [S15_vort_int.P2VMom2[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]
S15_vort_p2_v1_zif = [S15_vort_int.P2v1[n,S15.z_if_arg[n]] for n in range(0,S15.t_len)]

S20_vort_p2_s1_mean_zif_1 = [S20_vort_int_1.P2S1Mom1[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_vort_p2_s1_mean_zif_2 = [S20_vort_int_2.P2S1Mom1[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_vort_p2_s1_mean_zif = S20_vort_p2_s1_mean_zif_1 + S20_vort_p2_s1_mean_zif_2

S20_vort_p2_w_mean_zif_1 = [S20_vort_int_1.P2VMom1[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_vort_p2_w_mean_zif_2 = [S20_vort_int_2.P2VMom1[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_vort_p2_w_mean_zif = S20_vort_p2_w_mean_zif_1 + S20_vort_p2_w_mean_zif_2

S20_vort_p2_s1_var_zif_1 = [S20_vort_int_1.P2S1Mom2[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_vort_p2_s1_var_zif_2 = [S20_vort_int_2.P2S1Mom2[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_vort_p2_s1_var_zif = S20_vort_p2_s1_var_zif_1 + S20_vort_p2_s1_var_zif_2

S20_vort_p2_w_var_zif_1 = [S20_vort_int_1.P2VMom2[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_vort_p2_w_var_zif_2 = [S20_vort_int_2.P2VMom2[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_vort_p2_w_var_zif = S20_vort_p2_w_var_zif_1 + S20_vort_p2_w_var_zif_2

S20_vort_p2_v1_zif_1 = [S20_vort_int_1.P2v1[n,S20_1.z_if_arg[n]] for n in range(0,S20_1.t_len)]
S20_vort_p2_v1_zif_2 = [S20_vort_int_2.P2v1[n,S20_2.z_if_arg[n]] for n in range(0,S20_2.t_len)]
S20_vort_p2_v1_zif = S20_vort_p2_v1_zif_1 + S20_vort_p2_v1_zif_2

S25_vort_p2_s1_mean_zif = [S25_vort_int.P2S1Mom1[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]
S25_vort_p2_w_mean_zif = [S25_vort_int.P2VMom1[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]
S25_vort_p2_s1_var_zif = [S25_vort_int.P2S1Mom2[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]
S25_vort_p2_w_var_zif = [S25_vort_int.P2VMom2[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]
S25_vort_p2_v1_zif = [S25_vort_int.P2v1[n,S25.z_if_arg[n]] for n in range(0,S25.t_len)]

# time

time_1 = [S20_1.z_enc[n]/L0 for n in range(0,S20_1.t_len)]
time_2 = [S20_2.z_enc[n]/L0 for n in range(0,S20_2.t_len)]

time = time_1 + time_2

# scales
z_if = np.concatenate([S20_1.z_if,S20_2.z_if])
z_is = np.concatenate([S20_1.z_is,S20_2.z_if])
z_ig = np.concatenate([S20_1.z_ig,S20_2.z_ig])
z_enc = np.concatenate([S20_1.z_enc,S20_2.z_enc])
w_enc = (B0*z_enc)**(1./3.)

Delta_U = np.concatenate([S20_1.DeltaU,S20_2.DeltaU])

rho_bw_S0 = S0.Rsv/(S0.r2S*S0.Ryy)**0.5
rho_bw_turb_S0 = (S0_vort_int.P2v1[4:7,:]-S0_vort_int.P2S1Mom1[4:7,:]*S0_vort_int.P2VMom1[4:7,:])/(S0_vort_int.P2S1Mom2[4:7,:]*S0_vort_int.P2VMom2[4:7,:])**0.5
rho_bw_S05 = S05.Rsv/(S05.r2S*S05.Ryy)**0.5
rho_bw_turb_S05 = (S05_vort_int.P2v1[4:7,:]-S05_vort_int.P2S1Mom1[4:7,:]*S05_vort_int.P2VMom1[4:7,:])/(S05_vort_int.P2S1Mom2[4:7,:]*S05_vort_int.P2VMom2[4:7,:])**0.5
rho_bw_S10 = S10.Rsv/(S10.r2S*S10.Ryy)**0.5
rho_bw_turb_S10 = (S10_vort_int.P2v1[4:7,:]-S10_vort_int.P2S1Mom1[4:7,:]*S10_vort_int.P2VMom1[4:7,:])/(S10_vort_int.P2S1Mom2[4:7,:]*S10_vort_int.P2VMom2[4:7,:])**0.5
rho_bw_S15 = S15.Rsv/(S15.r2S*S15.Ryy)**0.5
rho_bw_turb_S15 = (S15_vort_int.P2v1[4:7,:]-S15_vort_int.P2S1Mom1[4:7,:]*S15_vort_int.P2VMom1[4:7,:])/(S15_vort_int.P2S1Mom2[4:7,:]*S15_vort_int.P2VMom2[4:7,:])**0.5
rho_bw_S20 = S20_1.Rsv/(S20_1.r2S*S20_1.Ryy)**0.5
rho_bw_S20_turb = (S20_vort_int_1.P2v1[4:,:]-S20_vort_int_1.P2S1Mom1[4:,:]*S20_vort_int_1.P2VMom1[4:,:])/(S20_vort_int_1.P2S1Mom2[4:,:]*S20_vort_int_1.P2VMom2[4:,:])**0.5
rho_bw_S25 = S25.Rsv/(S25.r2S*S25.Ryy)**0.5
rho_bw_turb_S25 = (S25_vort_int.P2v1[4:7,:]-S25_vort_int.P2S1Mom1[4:7,:]*S25_vort_int.P2VMom1[4:7,:])/(S25_vort_int.P2S1Mom2[4:7,:]*S25_vort_int.P2VMom2[4:7,:])**0.5
    
rho_bw_zif_S0 = np.array(S0_s1_flux_zif)/(np.array(S0_s1_var_zif)*np.array(S0_w_var_zif))**0.5
rho_bw_zif_S0_turb =  (np.array(S0_vort_p2_v1_zif)-np.array(S0_vort_p2_s1_mean_zif)*np.array(S0_vort_p2_w_mean_zif))/(np.array(S0_vort_p2_s1_var_zif)*np.array(S0_vort_p2_w_var_zif))**0.5
rho_bw_zif_S05 = np.array(S05_s1_flux_zif)/(np.array(S05_s1_var_zif)*np.array(S05_w_var_zif))**0.5
rho_bw_zif_S05_turb =  (np.array(S05_vort_p2_v1_zif)-np.array(S05_vort_p2_s1_mean_zif)*np.array(S05_vort_p2_w_mean_zif))/(np.array(S05_vort_p2_s1_var_zif)*np.array(S05_vort_p2_w_var_zif))**0.5
rho_bw_zif_S10 = np.array(S10_s1_flux_zif)/(np.array(S10_s1_var_zif)*np.array(S10_w_var_zif))**0.5
rho_bw_zif_S10_turb =  (np.array(S10_vort_p2_v1_zif)-np.array(S10_vort_p2_s1_mean_zif)*np.array(S10_vort_p2_w_mean_zif))/(np.array(S10_vort_p2_s1_var_zif)*np.array(S10_vort_p2_w_var_zif))**0.5
rho_bw_zif_S15 = np.array(S15_s1_flux_zif)/(np.array(S15_s1_var_zif)*np.array(S15_w_var_zif))**0.5
rho_bw_zif_S15_turb =  (np.array(S15_vort_p2_v1_zif)-np.array(S15_vort_p2_s1_mean_zif)*np.array(S15_vort_p2_w_mean_zif))/(np.array(S15_vort_p2_s1_var_zif)*np.array(S15_vort_p2_w_var_zif))**0.5
rho_bw_zif_S20 = np.array(S20_s1_flux_zif)/(np.array(S20_s1_var_zif)*np.array(S20_w_var_zif))**0.5
rho_bw_zif_S20_turb = (np.array(S20_vort_p2_v1_zif)-np.array(S20_vort_p2_s1_mean_zif)*np.array(S20_vort_p2_w_mean_zif))/(np.array(S20_vort_p2_s1_var_zif)*np.array(S20_vort_p2_w_var_zif))**0.5
rho_bw_zif_S25 = np.array(S25_s1_flux_zif)/(np.array(S25_s1_var_zif)*np.array(S25_w_var_zif))**0.5
rho_bw_zif_S25_turb =  (np.array(S25_vort_p2_v1_zif)-np.array(S25_vort_p2_s1_mean_zif)*np.array(S25_vort_p2_w_mean_zif))/(np.array(S25_vort_p2_s1_var_zif)*np.array(S25_vort_p2_w_var_zif))**0.5

#######################################################################
# Plot

blues = matplotlib.cm.get_cmap('Blues')
greens = matplotlib.cm.get_cmap('Greens')

f, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax3.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax3.tick_params(bottom=False,left=False)
ax1.set_xlim(-0.2,1.1)
ax1.set_ylim(-0.1,1.6)
ax2.set_ylim(0.3,1)
ax2.set_xlim(10,30)
ax1.plot(np.mean(S0_vort_int.int2[4:7,:],axis=0),S0.y/np.mean(S0.z_enc[4:7]),c=blues(0.5),label=r'$Fr_0=0$')
ax1.plot(np.mean(S05_vort_int.int2[4:7,:],axis=0),S05.y/np.mean(S05.z_enc[4:7]),c=blues(0.6),label=r'$Fr_0=5$')
ax1.plot(np.mean(S10_vort_int.int2[4:7,:],axis=0),S10.y/np.mean(S10.z_enc[4:7]),c=blues(0.7),label=r'$Fr_0=10$')
ax1.plot(np.mean(S15_vort_int.int2[4:7,:],axis=0),S15.y/np.mean(S15.z_enc[4:7]),c=blues(0.8),label=r'$Fr_0=15$')
ax1.plot(np.mean(S20_vort_int_1.int2[4:,:],axis=0),S20_1.y/np.mean(S20_1.z_enc[4:]),c=blues(0.9),label=r'$Fr_0=20$')
ax1.plot(np.mean(S25_vort_int.int2[4:7,:],axis=0),S25.y/np.mean(S25.z_enc[4:7]),c=blues(1.0),label=r'$Fr_0=25$')
ax2.plot(S0.z_enc[1:-1]/L0,runningmean(S0_vort_turbareafrac_zif,1),c=blues(0.5))
ax2.plot(S05.z_enc[1:-1]/L0,runningmean(S05_vort_turbareafrac_zif,1),c=blues(0.6))
ax2.plot(S10.z_enc[1:-1]/L0,runningmean(S10_vort_turbareafrac_zif,1),c=blues(0.7))
ax2.plot(S15.z_enc[1:-1]/L0,runningmean(S15_vort_turbareafrac_zif,1),c=blues(0.8))
ax2.plot(time[1:-1],runningmean(S20_vort_turbareafrac_zif,1),c=blues(0.9))
ax2.plot(S25.z_enc[1:-1]/L0,runningmean(S25_vort_turbareafrac_zif,1),c=blues(1.0))
ax3.plot(S05.DeltaU[1:-1]/(N*0.25*S05.z_enc[1:-1]),runningmean(np.array(S05_vort_turbareafrac_zif)/np.array(S0_vort_turbareafrac_zif),1),c=blues(0.6))
ax3.plot(S10.DeltaU[1:-1]/(N*0.25*S10.z_enc[1:-1]),runningmean(np.array(S10_vort_turbareafrac_zif)/np.array(S0_vort_turbareafrac_zif[8:]),1),c=blues(0.7))
ax3.plot(S15.DeltaU[1:-1]/(N*0.25*S15.z_enc[1:-1]),runningmean(np.array(S15_vort_turbareafrac_zif)/np.array(S0_vort_turbareafrac_zif[8:]),1),c=blues(0.8))
ax3.plot(Delta_U[1:-1]/(N*0.25*z_enc[1:-1]),runningmean(np.array(S20_vort_turbareafrac_zif)/np.array(S0_vort_turbareafrac_zif[8:]),1),c=blues(0.9))
ax3.plot(S25.DeltaU[1:-1]/(N*0.25*S25.z_enc[1:-1]),runningmean(np.array(S25_vort_turbareafrac_zif)/np.array(S0_vort_turbareafrac_zif[8:]),1),c=blues(1.0))
ax1.set_xlabel(r'$a_T$')
ax1.set_ylabel(r'$z/z_{enc}$')
ax1.set_title(r'(a)',fontsize=24,loc='left')
ax2.set_xlabel(r'$z_{enc}/L_0$')
ax2.set_ylabel(r'$(a_T)_{z_{i,f}}$')
ax2.set_title(r'(b)',fontsize=24,loc='left')
ax3.set_xlabel(r'$\Delta u/[N_0(\Delta z_i)_c]$')
ax3.set_ylabel(r'$(a_\mathrm{T})_{z_{i,f}}/[(a_\mathrm{T})_{z_{i,f}}]_c$')
ax3.set_title(r'(c)',fontsize=24,loc='left')
ax1.legend(loc='lower left',fontsize=18,handlelength=1,borderaxespad=0.1)
plt.tight_layout()
plt.savefig(opath+'turb_area_frac_extended.pdf',bbox_inches='tight')
plt.show()


f, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax3.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax3.tick_params(bottom=False,left=False)
ax1.set_ylim(0,1.4)
ax1.set_xlim(-0.3,1)
ax2.set_ylim(-0.05,0.3)
ax2.set_xlim(10,30)
ax1.plot(np.mean(rho_bw_turb_S0,axis=0),S0.y/np.mean(S0.z_enc[4:7]),c=blues(0.5))
ax1.plot(np.mean(rho_bw_turb_S05,axis=0),S05.y/np.mean(S05.z_enc[4:7]),c=blues(0.6))
ax1.plot(np.mean(rho_bw_turb_S10,axis=0),S10.y/np.mean(S10.z_enc[4:7]),c=blues(0.7))
ax1.plot(np.mean(rho_bw_turb_S15,axis=0),S15.y/np.mean(S15.z_enc[4:7]),c=blues(0.8))
ax1.plot(np.mean(rho_bw_S20_turb,axis=0),S20_1.y/np.mean(S20_1.z_enc[4:]),c=blues(0.9))
ax1.plot(np.mean(rho_bw_turb_S25,axis=0),S25.y/np.mean(S25.z_enc[4:7]),c=blues(1.0))
ax2.plot(S0.z_enc[1:-1]/L0,-runningmean(rho_bw_zif_S0_turb,1),c=blues(0.5),label=r'$Fr_0=0$')
ax2.plot(S05.z_enc[1:-1]/L0,-runningmean(rho_bw_zif_S05_turb,1),c=blues(0.6),label=r'$Fr_0=5$')
ax2.plot(S10.z_enc[1:-1]/L0,-runningmean(rho_bw_zif_S10_turb,1),c=blues(0.7),label=r'$Fr_0=10$')
ax2.plot(S15.z_enc[1:-1]/L0,-runningmean(rho_bw_zif_S15_turb,1),c=blues(0.8),label=r'$Fr_0=15$')
ax2.plot(time[1:-1],-runningmean(rho_bw_zif_S20_turb,1),c=blues(0.9),label=r'$Fr_0=20$')
ax2.plot(S25.z_enc[1:-1]/L0,-runningmean(rho_bw_zif_S25_turb,1),c=blues(1.0),label=r'$Fr_0=25$')
ax3.plot(S05.DeltaU[1:-1]/(N*0.25*S05.z_enc[1:-1]),runningmean(rho_bw_zif_S05_turb/rho_bw_zif_S0_turb,1),c=blues(0.6))
ax3.plot(S10.DeltaU[1:-1]/(N*0.25*S10.z_enc[1:-1]),runningmean(rho_bw_zif_S10_turb/rho_bw_zif_S0_turb[8:],1),c=blues(0.7))
ax3.plot(S15.DeltaU[1:-1]/(N*0.25*S15.z_enc[1:-1]),runningmean(rho_bw_zif_S15_turb/rho_bw_zif_S0_turb[8:],1),c=blues(0.8))
ax3.plot(Delta_U[1:-1]/(N*0.25*z_enc[1:-1]),runningmean(rho_bw_zif_S20_turb/rho_bw_zif_S0_turb[8:],1),c=blues(0.9))
ax3.plot(S25.DeltaU[1:-1]/(N*0.25*S25.z_enc[1:-1]),runningmean(rho_bw_zif_S25_turb/rho_bw_zif_S0_turb[8:],1),c=blues(1.0))
ax1.set_xlabel(r'$(\rho_{bw})_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_{enc}$')
ax1.set_title('(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_{enc}/L_0$')
ax2.set_ylabel(r'$-[(\rho_{bw})_\mathrm{T}]_{z_{i,f}}$')
ax2.set_title('(b)',fontsize=24,loc='left')
ax2.legend(loc='best',fontsize=18,handlelength=1,borderaxespad=0.1,ncol=2,columnspacing=1)
ax3.set_xlabel(r'$\Delta u/[N_0(\Delta z_i)_c]$')
ax3.set_ylabel(r'$[(\rho_{bw})_T]_{z_{i,f}}/[[(\rho_{bw})_T]_{z_{i,f}}]_c$')
ax3.set_title('(c)',fontsize=24,loc='left')
plt.tight_layout()
plt.savefig(opath+'rho_bw_extended.pdf',bbox_inches='tight')
plt.show()

f, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax3.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax3.tick_params(bottom=False,left=False)
ax1.set_xlim(-0.25,0.1)
ax1.set_ylim(1,1.4)
ax2.set_xlim(10,30)
ax2.set_ylim(0,0.25)
ax1.plot(np.mean(S0_vort_int.P2v1[4:7,:]-S0_vort_int.P2S1Mom1[4:7,:]*S0_vort_int.P2VMom1[4:7,:],axis=0)/B0,S0.y/np.mean(S0.z_enc[4:7]),c=blues(0.5))
ax1.plot(np.mean(S05_vort_int.P2v1[4:7,:]-S05_vort_int.P2S1Mom1[4:7,:]*S05_vort_int.P2VMom1[4:7,:],axis=0)/B0,S05.y/np.mean(S05.z_enc[4:7]),c=blues(0.6))
ax1.plot(np.mean(S10_vort_int.P2v1[4:7,:]-S10_vort_int.P2S1Mom1[4:7,:]*S10_vort_int.P2VMom1[4:7,:],axis=0)/B0,S10.y/np.mean(S10.z_enc[4:7]),c=blues(0.7))
ax1.plot(np.mean(S15_vort_int.P2v1[4:7,:]-S15_vort_int.P2S1Mom1[4:7,:]*S15_vort_int.P2VMom1[4:7,:],axis=0)/B0,S15.y/np.mean(S15.z_enc[4:7]),c=blues(0.8))
ax1.plot(np.mean(S20_vort_int_1.P2v1[4:,:]-S20_vort_int_1.P2S1Mom1[4:,:]*S20_vort_int_1.P2VMom1[4:,:],axis=0)/B0,S20_1.y/np.mean(S20_1.z_enc[4:]),c=blues(0.9))
ax1.plot(np.mean(S25_vort_int.P2v1[4:7,:]-S25_vort_int.P2S1Mom1[4:7,:]*S25_vort_int.P2VMom1[4:7,:],axis=0)/B0,S25.y/np.mean(S25.z_enc[4:7]),c=blues(1.0))
ax2.plot(S0.z_enc[1:-1]/L0,-runningmean(np.array(S0_vort_p2_v1_zif)-np.array(S0_vort_p2_s1_mean_zif)*np.array(S0_vort_p2_w_mean_zif),1)/B0,c=blues(0.5),label=r'$Fr_0=0$')
ax2.plot(S05.z_enc[1:-1]/L0,-runningmean(np.array(S05_vort_p2_v1_zif)-np.array(S05_vort_p2_s1_mean_zif)*np.array(S05_vort_p2_w_mean_zif),1)/B0,c=blues(0.6),label=r'$Fr_0=5$')
ax2.plot(S10.z_enc[1:-1]/L0,-runningmean(np.array(S10_vort_p2_v1_zif)-np.array(S10_vort_p2_s1_mean_zif)*np.array(S10_vort_p2_w_mean_zif),1)/B0,c=blues(0.7),label=r'$Fr_0=10$')
ax2.plot(S15.z_enc[1:-1]/L0,-runningmean(np.array(S15_vort_p2_v1_zif)-np.array(S15_vort_p2_s1_mean_zif)*np.array(S15_vort_p2_w_mean_zif),1)/B0,c=blues(0.8),label=r'$Fr_0=15$')
ax2.plot(time[1:-1],-runningmean(np.array(S20_vort_p2_v1_zif)-np.array(S20_vort_p2_s1_mean_zif)*np.array(S20_vort_p2_w_mean_zif),1)/B0,c=blues(0.9),label=r'$Fr_0=20$')
ax2.plot(S25.z_enc[1:-1]/L0,-runningmean(np.array(S25_vort_p2_v1_zif)-np.array(S25_vort_p2_s1_mean_zif)*np.array(S25_vort_p2_w_mean_zif),1)/B0,c=blues(1.0),label=r'$Fr_0=25$')
ax3.plot(S05.DeltaU[1:-1]/(N*0.25*S05.z_enc[1:-1]),runningmean((np.array(S05_vort_p2_v1_zif)-np.array(S05_vort_p2_s1_mean_zif)*np.array(S05_vort_p2_w_mean_zif))/(np.array(S0_vort_p2_v1_zif)-np.array(S0_vort_p2_s1_mean_zif)*np.array(S0_vort_p2_w_mean_zif)),1),c=blues(0.6))
ax3.plot(S10.DeltaU[1:-1]/(N*0.25*S10.z_enc[1:-1]),runningmean((np.array(S10_vort_p2_v1_zif)-np.array(S10_vort_p2_s1_mean_zif)*np.array(S10_vort_p2_w_mean_zif))/(np.array(S0_vort_p2_v1_zif[8:])-np.array(S0_vort_p2_s1_mean_zif[8:])*np.array(S0_vort_p2_w_mean_zif[8:])),1),c=blues(0.7))
ax3.plot(S15.DeltaU[1:-1]/(N*0.25*S15.z_enc[1:-1]),runningmean((np.array(S15_vort_p2_v1_zif)-np.array(S15_vort_p2_s1_mean_zif)*np.array(S15_vort_p2_w_mean_zif))/(np.array(S0_vort_p2_v1_zif[8:])-np.array(S0_vort_p2_s1_mean_zif[8:])*np.array(S0_vort_p2_w_mean_zif[8:])),1),c=blues(0.8))
ax3.plot(Delta_U[1:-1]/(N*0.25*z_enc[1:-1]),runningmean((np.array(S20_vort_p2_v1_zif)-np.array(S20_vort_p2_s1_mean_zif)*np.array(S20_vort_p2_w_mean_zif))/(np.array(S0_vort_p2_v1_zif[8:])-np.array(S0_vort_p2_s1_mean_zif[8:])*np.array(S0_vort_p2_w_mean_zif[8:])),1),c=blues(0.9))
ax3.plot(S25.DeltaU[1:-1]/(N*0.25*S25.z_enc[1:-1]),runningmean((np.array(S25_vort_p2_v1_zif)-np.array(S25_vort_p2_s1_mean_zif)*np.array(S25_vort_p2_w_mean_zif))/(np.array(S0_vort_p2_v1_zif[8:])-np.array(S0_vort_p2_s1_mean_zif[8:])*np.array(S0_vort_p2_w_mean_zif[8:])),1),c=blues(1.0))
ax1.set_xlabel(r'$\langle b^\prime w^\prime \rangle_\mathrm{T}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title('(a)',fontsize=20,loc='left')
ax2.set_xlabel(r'$z_{enc}/L_0$')
ax2.set_ylabel(r'$-(\langle b^\prime w^\prime \rangle_\mathrm{T})_{z_{i,f}}$')
ax2.set_title('(b)',fontsize=24,loc='left')
ax2.legend(loc='best',fontsize=18,handlelength=1,borderaxespad=0.1,ncol=2,columnspacing=1)
ax3.set_xlabel(r'$\Delta u/[N_0(\Delta z_i)_c]$')
ax3.set_ylabel(r'$(\langle b^\prime w^\prime\rangle_T)_{z_{i,f}}/[(\langle b^\prime w^\prime \rangle_T)_{z_{i,f}}]_c$')
ax3.set_title('(c)',fontsize=24,loc='left')
plt.tight_layout()
plt.savefig(opath+'s1_vflux_extended.pdf',bbox_inches='tight')
plt.show()

axes = plt.gca()
axes.grid(True)
plt.plot(S05.DeltaU[1:-1]/(N*0.25*S05.z_enc[1:-1]),runningmean(np.array(S05_vort_turbareafrac_zif)/np.array(S0_vort_turbareafrac_zif),1),c=blues(0.6),label=r'$Fr_0=5$')
plt.plot(S10.DeltaU[1:-1]/(N*0.25*S10.z_enc[1:-1]),runningmean(np.array(S10_vort_turbareafrac_zif)/np.array(S0_vort_turbareafrac_zif),1),c=blues(0.7),label=r'$Fr_0=10$')
plt.plot(S15.DeltaU[1:-1]/(N*0.25*S15.z_enc[1:-1]),runningmean(np.array(S15_vort_turbareafrac_zif)/np.array(S0_vort_turbareafrac_zif),1),c=blues(0.8),label=r'$Fr_0=15$')
plt.plot(Delta_U[1:-1]/(N*0.25*z_enc[1:-1]),runningmean(np.array(S20_vort_turbareafrac_zif)/np.array(S0_vort_turbareafrac_zif),1),c=blues(0.9),label=r'$Fr_0=20$')
plt.plot(S25.DeltaU[1:-1]/(N*0.25*S25.z_enc[1:-1]),runningmean(np.array(S25_vort_turbareafrac_zif)/np.array(S0_vort_turbareafrac_zif),1),c=blues(1.0),label=r'$Fr_0=25$')
plt.xlabel(r'$\Delta u/[N_0(\Delta z_i)_c]$')
plt.ylabel(r'$(a_\mathrm{T})_{z_{i,f}}/((a_\mathrm{T})_{z_{i,f}})_c$')
plt.legend(loc='best',handlelength=1,borderaxespad=0.1)
plt.tight_layout()
plt.savefig(opath+'turb_area_frac_deltaU.pdf',bbox_inches='tight')
plt.show()

axes = plt.gca()
axes.grid(True)
plt.ylim(0,1)
plt.plot(S05.DeltaU[1:-1]/(N*0.25*S05.z_enc[1:-1]),runningmean(rho_bw_zif_S05_turb/rho_bw_zif_S0_turb,1),c=blues(0.6),label=r'$Fr_0=5$')
plt.plot(S10.DeltaU[1:-1]/(N*0.25*S10.z_enc[1:-1]),runningmean(rho_bw_zif_S10_turb/rho_bw_zif_S0_turb,1),c=blues(0.7),label=r'$Fr_0=10$')
plt.plot(S15.DeltaU[1:-1]/(N*0.25*S15.z_enc[1:-1]),runningmean(rho_bw_zif_S15_turb/rho_bw_zif_S0_turb,1),c=blues(0.8),label=r'$Fr_0=15$')
plt.plot(Delta_U[1:-1]/(N*0.25*z_enc[1:-1]),runningmean(rho_bw_zif_S20_turb/rho_bw_zif_S0_turb,1),c=blues(0.9),label=r'$Fr_0=20$')
plt.plot(S25.DeltaU[1:-1]/(N*0.25*S25.z_enc[1:-1]),runningmean(rho_bw_zif_S25_turb/rho_bw_zif_S0_turb,1),c=blues(1.0),label=r'$Fr_0=25$')
plt.xlabel(r'$\Delta u/[N_0(\Delta z_i)_c]$')
plt.ylabel(r'$(\rho_{bw})_{z_{i,f}}/((\rho_{bw})_{z_{i,f}})_c$')
plt.legend(loc='best',fontsize=18,handlelength=1,borderaxespad=0.1,ncol=2,columnspacing=1)
plt.tight_layout()
plt.savefig(opath+'rho_bw_deltaU.pdf',bbox_inches='tight')
plt.show()


axes = plt.gca()
axes.grid(True)
plt.ylim(0,1)
plt.plot(S10.DeltaU[1:-1]/(N*0.25*S10.z_enc[1:-1]),runningmean((np.array(S10_vort_p2_v1_zif)-np.array(S10_vort_p2_s1_mean_zif)*np.array(S10_vort_p2_w_mean_zif))/(np.array(S0_vort_p2_v1_zif)-np.array(S0_vort_p2_s1_mean_zif)*np.array(S0_vort_p2_w_mean_zif)),1),c=blues(0.7),label=r'$Fr_0=10$')
plt.plot(S15.DeltaU[1:-1]/(N*0.25*S15.z_enc[1:-1]),runningmean((np.array(S15_vort_p2_v1_zif)-np.array(S15_vort_p2_s1_mean_zif)*np.array(S15_vort_p2_w_mean_zif))/(np.array(S0_vort_p2_v1_zif)-np.array(S0_vort_p2_s1_mean_zif)*np.array(S0_vort_p2_w_mean_zif)),1),c=blues(0.8),label=r'$Fr_0=15$')
plt.plot(Delta_U[1:-1]/(N*0.25*z_enc[1:-1]),runningmean((np.array(S20_vort_p2_v1_zif)-np.array(S20_vort_p2_s1_mean_zif)*np.array(S20_vort_p2_w_mean_zif))/(np.array(S0_vort_p2_v1_zif)-np.array(S0_vort_p2_s1_mean_zif)*np.array(S0_vort_p2_w_mean_zif)),1),c=blues(0.9),label=r'$Fr_0=20$')
plt.xlabel(r'$\Delta u/[N_0(\Delta z_i)_c]$')
plt.ylabel(r'$(\langle b^\prime w^\prime \rangle_\mathrm{T})_{z_{i,f}}/((\langle b^\prime w^\prime \rangle_\mathrm{T})_{z_{i,f}})_c$')
plt.legend(loc='best',handlelength=1)
plt.tight_layout()
plt.savefig(opath+'s1_vflux_deltaU.pdf',bbox_inches='tight')
plt.show()

axes = plt.gca()
axes.grid(True)
plt.ylim(0,2)
plt.plot(S10.DeltaU[1:-1]/(N*0.25*S10.z_enc[1:-1]),runningmean((np.array(S10_vort_turbareafrac_zif)*(np.array(S10_vort_p2_v1_zif)-np.array(S10_vort_p2_s1_mean_zif)*np.array(S10_vort_p2_w_mean_zif)))/(np.array(S0_vort_turbareafrac_zif)*(np.array(S0_vort_p2_v1_zif)-np.array(S0_vort_p2_s1_mean_zif)*np.array(S0_vort_p2_w_mean_zif))),1),c=blues(0.7),label=r'$Fr_0=10$')
plt.plot(S15.DeltaU[1:-1]/(N*0.25*S15.z_enc[1:-1]),runningmean((np.array(S15_vort_turbareafrac_zif)*(np.array(S15_vort_p2_v1_zif)-np.array(S15_vort_p2_s1_mean_zif)*np.array(S15_vort_p2_w_mean_zif)))/(np.array(S0_vort_turbareafrac_zif)*(np.array(S0_vort_p2_v1_zif)-np.array(S0_vort_p2_s1_mean_zif)*np.array(S0_vort_p2_w_mean_zif))),1),c=blues(0.8),label=r'$Fr_0=15$')
plt.plot(Delta_U[1:-1]/(N*0.25*z_enc[1:-1]),runningmean((np.array(S20_vort_turbareafrac_zif)*(np.array(S20_vort_p2_v1_zif)-np.array(S20_vort_p2_s1_mean_zif)*np.array(S20_vort_p2_w_mean_zif)))/(np.array(S0_vort_turbareafrac_zif)*(np.array(S0_vort_p2_v1_zif)-np.array(S0_vort_p2_s1_mean_zif)*np.array(S0_vort_p2_w_mean_zif))),1),c=blues(0.9),label=r'$Fr_0=20$')
plt.xlabel(r'$\Delta u/[N_0(\Delta z_i)_c]$')
plt.ylabel(r'$(a_\mathrm{T}\langle b^\prime w^\prime\rangle )_{z_{i,f}}/((a_\mathrm{T}\langle b^\prime w^\prime \rangle_\mathrm{T})_{z_{i,f}})_c$')
plt.legend(loc='best',handlelength=1)
plt.tight_layout()
plt.savefig(opath+'aT_s1_vflux_deltaU.pdf',bbox_inches='tight')
plt.show()

axes = plt.gca()
axes.grid(True)
#axes.set_xlim(0,3)
axes.set_ylim(0.9,1.5)
plt.plot(S10.DeltaU[1:-1]/(N*0.25*S10.z_enc[1:-1]),runningmean(np.array(S10.z_if)/np.array(S10.z_enc),1),c=blues(0.7),label=r'$Fr_0=10$')
plt.plot(S15.DeltaU[1:-1]/(N*0.25*S15.z_enc[1:-1]),runningmean(np.array(S15.z_if)/np.array(S15.z_enc),1),c=blues(0.8),label=r'$Fr_0=15$')
plt.plot(Delta_U[1:-1]/(N*0.25*z_enc[1:-1]),runningmean(z_if/z_enc,1),c=blues(0.9),label=r'$Fr_0=20$')
# plt.plot(S10.DeltaU[1:-1]/(N*0.25*S10.z_enc[1:-1]),runningmean(np.array(S10.z_is)/np.array(S10.z_enc),1),c=blues(0.7))
# plt.plot(S15.DeltaU[1:-1]/(N*0.25*S15.z_enc[1:-1]),runningmean(np.array(S15.z_is)/np.array(S15.z_enc),1),c=blues(0.8))
# plt.plot(Delta_U[1:-1]/(N*0.25*z_enc[1:-1]),runningmean(z_is/z_enc,1),c=blues(0.9))
plt.xlabel(r'$\Delta u/[N_0(\Delta z_i)_c]$')
plt.ylabel(r'$z_{i,f}/z_\mathrm{enc}$')
plt.legend(loc='best',handlelength=1)
plt.tight_layout()
plt.savefig(opath+'z_if_deltaU.pdf',bbox_inches='tight')
plt.show()

axes = plt.gca()
axes.grid(True)
axes.set_ylim(0,1.6)
plt.plot(np.sqrt(np.mean(S0.r2S[4:7,:],axis=0))/(B0/(B0*np.mean(S0.z_enc[4:7]))**(1./3.)),S0.y/np.mean(S0.z_enc[4:7]),c=blues(0.5),label=r'$Fr_0=0$')
plt.plot(np.sqrt(np.mean(S20_1.r2S[4:7,:],axis=0))/(B0/(np.mean(w_enc[4:7]))),S20_2.y/np.mean(z_enc[4:7]),c=blues(0.9),label=r'$Fr_0=20$')
plt.xlabel(r'$b_\mathrm{rms}/b_*$')
plt.ylabel(r'$z/z_\mathrm{enc}$')
plt.legend(loc='best',handlelength=1)
plt.tight_layout()
plt.savefig(opath+'b_rms.pdf',bbox_inches='tight')
plt.show()



f, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax1.set_xlim(1,1.4)
ax2.set_xlim(-0.15,0.15)
ax1.set_ylim(1,1.4)
ax2.set_ylim(1,1.4)
ax1.plot(np.mean(NS42_vort_int_1.P2S1Mom1[-4:-1,:],axis=0)/np.mean(N**2*NS42_1.z_enc[-4:-1]),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),label=r'$Fr_0=0$')
ax1.plot(np.mean(S20_vort_int.P2S1Mom1[4:7,:],axis=0)/np.mean(N**2*S20.z_enc[4:7]),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),label=r'$Fr_0=20$')
ax1.plot(np.mean(NS42_vort_int_1.P1S1Mom1[-4:-1,:],axis=0)/np.mean(N**2*NS42_1.z_enc[-4:-1]),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),ls='--')
ax1.plot(np.mean(S20_vort_int.P1S1Mom1[4:7,:],axis=0)/np.mean(N**2*S20.z_enc[4:7]),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),ls='--')
ax1.plot(S20.y/np.mean(S20.z_enc[4:7]),S20.y/np.mean(S20.z_enc[4:7]),'k:')
ax2.plot(np.mean(NS42_vort_int_1.P2VMom1[-4:-1,:],axis=0)/(np.mean(NS42_1.z_enc[-4:-1])*B0)**(1./3.),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5))
ax2.plot(np.mean(S20_vort_int.P2VMom1[4:7,:],axis=0)/(np.mean(S20.z_enc[4:7])*B0)**(1./3.),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9))
ax2.plot(np.mean(NS42_vort_int_1.P1VMom1[-4:-1,:],axis=0)/(np.mean(NS42_1.z_enc[-4:-1])*B0)**(1./3.),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),ls='--')
ax2.plot(np.mean(S20_vort_int.P1VMom1[4:7,:],axis=0)/(np.mean(S20.z_enc[4:7])*B0)**(1./3.),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9),ls='--')
black_solid = mlines.Line2D([],[],c='k',label='Turbulent')
black_dashed = mlines.Line2D([],[],c='k',ls='--',label='Non-turbulent')
black_dot = mlines.Line2D([],[],c='k',ls=':',label=r'$N_0^2z$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax2.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax2.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$\langle b \rangle / b_\mathrm{enc}$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title(r'(a)',fontsize=24,loc='left')
ax2.set_xlabel(r'$\langle w \rangle / w_*$')
ax2.set_ylabel(r'$z/z_\mathrm{enc}$')
ax2.set_title(r'(b)',fontsize=24,loc='left')
leg1 = ax1.legend(loc='upper left',fontsize=18,borderaxespad=0.1,handlelength=1)
leg2 = ax1.legend(handles=[black_solid,black_dashed,black_dot],loc='lower right',fontsize=18,borderaxespad=0.1,handlelength=1,labelspacing=0.3)
ax1.add_artist(leg1)
plt.tight_layout()
plt.savefig(opath+'b_mean_w_mean.pdf',bbox_inches='tight')
plt.show()


f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(10,10))
ax1.grid(True,linewidth=1.5)
ax2.grid(True,linewidth=1.5)
ax3.grid(True,linewidth=1.5)
ax4.grid(True,linewidth=1.5)
ax1.tick_params(bottom=False,left=False)
ax2.tick_params(bottom=False,left=False)
ax3.tick_params(bottom=False,left=False)
ax4.tick_params(bottom=False,left=False)
ax1.set_xlim(0,2.5)
ax1.set_ylim(1,1.4)
ax2.set_xlim(15,30)
ax2.set_ylim(0,2)
ax3.set_ylim(1,1.4)
ax3.set_xlim(0,0.6)
ax4.set_ylim(0.2,0.6)
ax4.set_xlim(15,30)
ax1.plot(np.mean(np.sqrt(NS42_vort_int_1.P2S1Mom2[-4:-1,:]),axis=0)/(N**2*L0),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5),label=r'$Fr_0=0$')
ax1.plot(np.mean(np.sqrt(S20_vort_int.P2S1Mom2[4:7,:]),axis=0)/(N**2*L0),S20.y/np.mean(S20.z_enc[4:7],axis=0),c=blues(0.9),label=r'$Fr_0=20$')
ax2.plot(time[1:-1],runningmean(np.sqrt(NS42_vort_p2_s1_var_zif),1)/(N**2*L0),c=blues(0.5))
ax2.plot(S20.z_enc[1:-1]/L0,runningmean(np.sqrt(S20_vort_p2_s1_var_zif),1)/(N**2*L0),c=blues(0.9))
ax3.plot(np.mean(np.sqrt(NS42_vort_int_1.P2VMom2[-4:-1,:]),axis=0)/np.mean(w_enc[4:7]),NS42_1.y/np.mean(NS42_1.z_enc[-4:-1]),c=blues(0.5))
ax3.plot(np.mean(np.sqrt(S20_vort_int.P2VMom2[4:7,:]),axis=0)/np.mean((B0*S20.z_enc[4:7])**(1./3.)),S20.y/np.mean(S20.z_enc[4:7]),c=blues(0.9))
ax4.plot(time[1:-1],runningmean(np.sqrt(NS42_vort_p2_w_var_zif)/w_enc,1),c=blues(0.5),label=r'$Fr_0=0$')
ax4.plot(S20.z_enc[1:-1]/L0,runningmean(np.sqrt(S20_vort_p2_w_var_zif)/(B0*S20.z_enc)**(1./3.),1),c=blues(0.9),label=r'$Fr_0=20$')
ax1.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax1.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax3.axhline(np.mean(NS42_1.z_if[-4:-1]/NS42_1.z_enc[-4:-1]),0,0.05,c=blues(0.5))
ax3.axhline(np.mean(S20.z_if[4:7]/S20.z_enc[4:7]),0,0.05,c=blues(0.9))
ax1.set_xlabel(r'$(b_\mathrm{rms})_\mathrm{T} /(N^2L_0)$')
ax1.set_ylabel(r'$z/z_\mathrm{enc}$')
ax1.set_title('(a)',fontsize=24,loc='left')
ax2.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax2.set_ylabel(r'$((b_\mathrm{rms})_\mathrm{T})_{z_{i,f}}/(N^2L_0)$')
ax2.set_title('(b)',fontsize=24,loc='left')
ax3.set_xlabel(r'$(w_\mathrm{rms})_\mathrm{T} / w_*$')
ax3.set_ylabel(r'$z/z_\mathrm{enc}$')
ax3.set_title('(c)',fontsize=24,loc='left')
ax4.set_xlabel(r'$z_\mathrm{enc}/L_0$')
ax4.set_ylabel(r'$((w_\mathrm{rms})_\mathrm{T})_{z_{i,f}}/w_*$')
ax4.set_title('(d)',fontsize=24,loc='left')
ax1.legend(loc='best',fontsize=24,handlelength=1,borderaxespad=0.1)
plt.tight_layout()
plt.savefig(opath+'b_rms_w_rms.pdf',bbox_inches='tight')
plt.show()


axes = plt.gca()
axes.grid(True)
plt.xlim(-0.2,0.05)
plt.ylim(1,1.4)
plt.plot(np.mean(S0_vort_int.int2[4:7,:]*(S0_vort_int.P2v1[4:7,:]-S0_vort_int.P2S1Mom1[4:7,:]*S0_vort_int.P2VMom1[4:7,:])/B0,axis=0),S0.y/np.mean(S0.z_enc[4:7]),c=blues(0.5))
plt.plot(np.mean(S05_vort_int.int2[4:7,:]*(S05_vort_int.P2v1[4:7,:]-S05_vort_int.P2S1Mom1[4:7,:]*S05_vort_int.P2VMom1[4:7,:])/B0,axis=0),S05.y/np.mean(S05.z_enc[4:7]),c=blues(0.6))
plt.plot(np.mean(S10_vort_int.int2[4:7,:]*(S10_vort_int.P2v1[4:7,:]-S10_vort_int.P2S1Mom1[4:7,:]*S10_vort_int.P2VMom1[4:7,:])/B0,axis=0),S10.y/np.mean(S10.z_enc[4:7]),c=blues(0.7))
plt.plot(np.mean(S15_vort_int.int2[4:7,:]*(S15_vort_int.P2v1[4:7,:]-S15_vort_int.P2S1Mom1[4:7,:]*S15_vort_int.P2VMom1[4:7,:])/B0,axis=0),S15.y/np.mean(S15.z_enc[4:7]),c=blues(0.8))
plt.plot(np.mean(S20_vort_int_1.int2[4:,:]*(S20_vort_int_1.P2v1[4:,:]-S20_vort_int_1.P2S1Mom1[4:,:]*S20_vort_int_1.P2VMom1[4:,:])/B0,axis=0),S20_1.y/np.mean(S20_1.z_enc[4:]),c=blues(0.9))
plt.plot(np.mean(S25_vort_int.int2[4:7,:]*(S25_vort_int.P2v1[4:7,:]-S25_vort_int.P2S1Mom1[4:7,:]*S25_vort_int.P2VMom1[4:7,:])/B0,axis=0),S25.y/np.mean(S25.z_enc[4:7]),c=blues(1.0))
plt.plot(np.mean(S0_vort_int.int1[4:7,:]*(S0_vort_int.P1v1[4:7,:]-S0_vort_int.P1S1Mom1[4:7,:]*S0_vort_int.P1VMom1[4:7,:])/B0,axis=0),S0.y/np.mean(S0.z_enc[4:7]),c=blues(0.5),ls='--')
plt.plot(np.mean(S05_vort_int.int1[4:7,:]*(S05_vort_int.P1v1[4:7,:]-S05_vort_int.P1S1Mom1[4:7,:]*S05_vort_int.P1VMom1[4:7,:])/B0,axis=0),S05.y/np.mean(S05.z_enc[4:7]),c=blues(0.6),ls='--')
plt.plot(np.mean(S10_vort_int.int1[4:7,:]*(S10_vort_int.P1v1[4:7,:]-S10_vort_int.P1S1Mom1[4:7,:]*S10_vort_int.P1VMom1[4:7,:])/B0,axis=0),S10.y/np.mean(S10.z_enc[4:7]),c=blues(0.7),ls='--')
plt.plot(np.mean(S15_vort_int.int1[4:7,:]*(S15_vort_int.P1v1[4:7,:]-S15_vort_int.P1S1Mom1[4:7,:]*S15_vort_int.P1VMom1[4:7,:])/B0,axis=0),S15.y/np.mean(S15.z_enc[4:7]),c=blues(0.8),ls='--')
plt.plot(np.mean(S20_vort_int_1.int1[4:,:]*(S20_vort_int_1.P1v1[4:,:]-S20_vort_int_1.P1S1Mom1[4:,:]*S20_vort_int_1.P1VMom1[4:,:])/B0,axis=0),S20_1.y/np.mean(S20_1.z_enc[4:]),c=blues(0.9),ls='--')
plt.plot(np.mean(S25_vort_int.int1[4:7,:]*(S25_vort_int.P1v1[4:7,:]-S25_vort_int.P1S1Mom1[4:7,:]*S25_vort_int.P1VMom1[4:7,:])/B0,axis=0),S25.y/np.mean(S25.z_enc[4:7]),c=blues(1.0),ls='--')
plt.plot(np.mean(S0_vort_int.int2[4:7,:]*S0_vort_int.int1[4:7,:]*(S0_vort_int.P2S1Mom1[4:7,:]-S0_vort_int.P1S1Mom1[4:7,:])*(S0_vort_int.P2VMom1[4:7,:]-S0_vort_int.P1VMom1[4:7,:])/B0,axis=0),S0.y/np.mean(S0.z_enc[4:7]),c=blues(0.5),ls='-.')
plt.plot(np.mean(S05_vort_int.int2[4:7,:]*S05_vort_int.int1[4:7,:]*(S05_vort_int.P2S1Mom1[4:7,:]-S05_vort_int.P1S1Mom1[4:7,:])*(S05_vort_int.P2VMom1[4:7,:]-S05_vort_int.P1VMom1[4:7,:])/B0,axis=0),S05.y/np.mean(S05.z_enc[4:7]),c=blues(0.6),ls='-.')
plt.plot(np.mean(S10_vort_int.int2[4:7,:]*S10_vort_int.int1[4:7,:]*(S10_vort_int.P2S1Mom1[4:7,:]-S10_vort_int.P1S1Mom1[4:7,:])*(S10_vort_int.P2VMom1[4:7,:]-S10_vort_int.P1VMom1[4:7,:])/B0,axis=0),S10.y/np.mean(S10.z_enc[4:7]),c=blues(0.7),ls='-.')
plt.plot(np.mean(S15_vort_int.int2[4:7,:]*S15_vort_int.int1[4:7,:]*(S15_vort_int.P2S1Mom1[4:7,:]-S15_vort_int.P1S1Mom1[4:7,:])*(S15_vort_int.P2VMom1[4:7,:]-S15_vort_int.P1VMom1[4:7,:])/B0,axis=0),S15.y/np.mean(S15.z_enc[4:7]),c=blues(0.8),ls='-.')
plt.plot(np.mean(S20_vort_int_1.int2[4:,:]*S20_vort_int_1.int1[4:,:]*(S20_vort_int_1.P2S1Mom1[4:,:]-S20_vort_int_1.P1S1Mom1[4:,:])*(S20_vort_int_1.P2VMom1[4:,:]-S20_vort_int_1.P1VMom1[4:,:])/B0,axis=0),S20_1.y/np.mean(S20_1.z_enc[4:]),c=blues(0.9),ls='-.')
plt.plot(np.mean(S25_vort_int.int2[4:7,:]*S25_vort_int.int1[4:7,:]*(S25_vort_int.P2S1Mom1[4:7,:]-S25_vort_int.P1S1Mom1[4:7,:])*(S25_vort_int.P2VMom1[4:7,:]-S25_vort_int.P1VMom1[4:7,:])/B0,axis=0),S25.y/np.mean(S25.z_enc[4:7]),c=blues(0.5),ls='-.')
plt.plot(np.mean(S0.Rsv[4:7,:]/B0,axis=0),S0.y/np.mean(S0.z_enc[4:7]),c=blues(0.5),ls=':')
plt.plot(np.mean(S05.Rsv[4:7,:]/B0,axis=0),S05.y/np.mean(S05.z_enc[4:7]),c=blues(0.6),ls=':')
plt.plot(np.mean(S10.Rsv[4:7,:]/B0,axis=0),S10.y/np.mean(S10.z_enc[4:7]),c=blues(0.7),ls=':')
plt.plot(np.mean(S15.Rsv[4:7,:]/B0,axis=0),S15.y/np.mean(S15.z_enc[4:7]),c=blues(0.8),ls=':')
plt.plot(np.mean(S20_1.Rsv[4:,:]/B0,axis=0),S20_1.y/np.mean(S20_1.z_enc[4:]),c=blues(0.9),ls=':')
plt.plot(np.mean(S25.Rsv[4:7,:]/B0,axis=0),S25.y/np.mean(S25.z_enc[4:7]),c=blues(1.0),ls=':')
# plt.xlabel(r'$f /B_0$')
# plt.ylabel(r'$z/z_{enc}$')
# black_dot = mlines.Line2D([],[],c='k',ls=':',label=r'$\langle b^\prime w^\prime \rangle$')
# black_solid = mlines.Line2D([],[],c='k',label=r'$a_T\langle b^\prime w^\prime \rangle_T$')
# black_dashed = mlines.Line2D([],[],c='k',ls='--',label=r'$a_{NT}\langle b^\prime w^\prime \rangle_{NT}$')
# black_dashdot = mlines.Line2D([],[],c='k',ls='-.',label=r'$a_{T}a_{NT}(\langle b \rangle_{T} - \langle b \rangle_{NT})(\langle w \rangle_{T} - \langle w \rangle_{NT})$')
# plt.legend(handles=[black_dot,black_solid,black_dashed,black_dashdot],loc='upper right',fontsize=20,frameon=False,handlelength=1)
plt.tight_layout()
plt.savefig(opath+'s1_vflux_contributions.pdf',bbox_inches='tight')
plt.show()



axes = plt.gca()
axes.grid(True)
plt.xlim(0,1)
plt.ylim(1,1.4)
#plt.plot(np.mean(S0_vort_int.int2[4:7,:]*(S0_vort_int.P2v1[4:7,:]-S0_vort_int.P2S1Mom1[4:7,:]*S0_vort_int.P2VMom1[4:7,:])/S0.Rsv[4:7,:],axis=0),S0.y/np.mean(S0.z_enc[4:7]),c=blues(0.5))
plt.plot(np.mean(S05_vort_int.int2[4:7,:]*(S05_vort_int.P2v1[4:7,:]-S05_vort_int.P2S1Mom1[4:7,:]*S05_vort_int.P2VMom1[4:7,:])/S05.Rsv[4:7,:],axis=0),S05.y/np.mean(S05.z_enc[4:7]),c=blues(0.6))
plt.plot(np.mean(S10_vort_int.int2[4:7,:]*(S10_vort_int.P2v1[4:7,:]-S10_vort_int.P2S1Mom1[4:7,:]*S10_vort_int.P2VMom1[4:7,:])/S10.Rsv[4:7,:],axis=0),S10.y/np.mean(S10.z_enc[4:7]),c=blues(0.7))
#plt.plot(np.mean(S15_vort_int.int2[4:7,:]*(S15_vort_int.P2v1[4:7,:]-S15_vort_int.P2S1Mom1[4:7,:]*S15_vort_int.P2VMom1[4:7,:])/S15.Rsv[4:7,:],axis=0),S15.y/np.mean(S15.z_enc[4:7]),c=blues(0.8))
#plt.plot(np.mean(S20_vort_int_1.int2[4:,:]*(S20_vort_int_1.P2v1[4:,:]-S20_vort_int_1.P2S1Mom1[4:,:]*S20_vort_int_1.P2VMom1[4:,:])/S20_1.Rsv[4:,:],axis=0),S20_1.y/np.mean(S20_1.z_enc[4:]),c=blues(0.9))
#plt.plot(np.mean(S25_vort_int.int2[4:7,:]*(S25_vort_int.P2v1[4:7,:]-S25_vort_int.P2S1Mom1[4:7,:]*S25_vort_int.P2VMom1[4:7,:])/S25.Rsv[4:7,:],axis=0),S25.y/np.mean(S25.z_enc[4:7]),c=blues(1.0))
plt.plot(np.mean(S05_vort_int.int1[4:7,:]*(S05_vort_int.P1v1[4:7,:]-S05_vort_int.P1S1Mom1[4:7,:]*S05_vort_int.P1VMom1[4:7,:])/S05.Rsv[4:7,:],axis=0),S05.y/np.mean(S05.z_enc[4:7]),c=blues(0.6),ls='--')
plt.plot(np.mean(S10_vort_int.int1[4:7,:]*(S10_vort_int.P1v1[4:7,:]-S10_vort_int.P1S1Mom1[4:7,:]*S10_vort_int.P1VMom1[4:7,:])/S10.Rsv[4:7,:],axis=0),S10.y/np.mean(S10.z_enc[4:7]),c=blues(0.7),ls='--')
plt.plot(np.mean(S05_vort_int.int2[4:7,:]*S05_vort_int.int1[4:7,:]*(S05_vort_int.P2S1Mom1[4:7,:]-S05_vort_int.P1S1Mom1[4:7,:])*(S05_vort_int.P2VMom1[4:7,:]-S05_vort_int.P1VMom1[4:7,:])/S05.Rsv[4:7,:],axis=0),S05.y/np.mean(S05.z_enc[4:7]),c=blues(0.6),ls='-.')
plt.plot(np.mean(S10_vort_int.int2[4:7,:]*S10_vort_int.int1[4:7,:]*(S10_vort_int.P2S1Mom1[4:7,:]-S10_vort_int.P1S1Mom1[4:7,:])*(S10_vort_int.P2VMom1[4:7,:]-S10_vort_int.P1VMom1[4:7,:])/S10.Rsv[4:7,:],axis=0),S10.y/np.mean(S10.z_enc[4:7]),c=blues(0.7),ls='-.')
plt.xlabel(r'$f / \langle b^\prime w^\prime \rangle$')
plt.ylabel(r'$z/z_{enc}$')
plt.tight_layout()
plt.savefig(opath+'s1_vflux_contributions_percent.pdf',bbox_inches='tight')
plt.show()


axes = plt.gca()
axes.grid(True)
plt.plot(S0.z_enc[1:-1]/L0,-runningmean(S0_s1_flux_zif,1)/B0,c=blues(0.5))
plt.plot(S05.z_enc[1:-1]/L0,-runningmean(S05_s1_flux_zif,1)/B0,c=blues(0.6))
plt.plot(S10.z_enc[1:-1]/L0,-runningmean(S10_s1_flux_zif,1)/B0,c=blues(0.7))
plt.plot(S15.z_enc[1:-1]/L0,-runningmean(S15_s1_flux_zif,1)/B0,c=blues(0.8))
plt.plot(time[1:-1],-runningmean(S20_s1_flux_zif,1)/B0,c=blues(0.9))
plt.plot(S25.z_enc[1:-1]/L0,-runningmean(S25_s1_flux_zif,1)/B0,c=blues(1.0))
plt.xlabel(r'$z_{enc}/L_0$')
plt.ylabel(r'$-(\langle b^\prime w^\prime \rangle)_{z_{i,f}}/B_0$')
plt.tight_layout()
plt.savefig(opath+'s1_vflux_total.pdf',bbox_inches='tight')
plt.show()
