from pathlib import Path 
import xmitgcm 
from matplotlib import pyplot as plt 
import numpy as np 
import xarray as xr 
from xgcm import Grid 
from scipy.ndimage import uniform_filter 
import matplotlib.animation as animation 
from IPython.display import HTML 
import pandas as pd 
import matplotlib.dates as mdates 
import matplotlib.lines as mlines 
import netCDF4 
# import zarr

import time 
import funcs.plotting as plot 
import funcs.velocity as decomp 
import funcs.IterativeSolvers as IS 
import funcs.analysis as an 
import warnings 
warnings.filterwarnings('ignore') 
%load_ext autoreload 
%autoreload 2

# data_dir = Path('/storage/mloch/arctic36km/run_final')
data_dir = Path('/storage/mloch//narvaltest/run_narvalJuly2026')

# Daily
dsSI   = xmitgcm.open_mdsdataset(data_dir,ref_date='1992-1-1 0:0:0', delta_t=1800, swap_dims=False, prefix = ['ds_SI']) #.isel(time=slice(50,None)) 
dsvert  = xmitgcm.open_mdsdataset(data_dir,ref_date='1992-1-1 0:0:0', delta_t=1800, swap_dims=False, prefix = ['ds_vert'])# .isel(time=slice(50,None))
dsATM  = xmitgcm.open_mdsdataset(data_dir,ref_date='1992-1-1 0:0:0', delta_t=1800, swap_dims=False, prefix = ['ds_atm'])# .isel(time=slice(50,None))

grid = Grid(dsvert,['X','Y'])
dsSI

metrics = {
    ('X',): ['dxC', 'dxG'],
    ('Y',): ['dyC', 'dyG'],
    ('Z',): ['drC', 'drF'],
    ('X', 'Y'): ['rA', 'rAw', 'rAs', 'rAz']}

grid = Grid(dsSI, metrics=metrics)


## DEFINE FOR EACH FUN 
# Ex to get the mean of the March means: group = 'time.month' and time =3 
config = 'arctic36km'
sim = 'converged'

group = 'time.month'
time = 2

if time ==0:
    month = 'Jan_'
if time ==1:
    month = 'Feb_'
if time ==2:
    month = 'Mars_'
if time ==3:
    month = 'Apr_'
if time ==4:
    month = 'May_'
if time ==5:
    month = 'Jun_'
if time ==6:
    month = 'Jul_'
if time ==7:
    month = 'Aug_'
if time ==8:
    month = 'Sept_'
if time ==9:
    month = 'Oct_'
if time ==10:
    month = 'Nov_'
if time ==11:
    month = 'Dec_'

    # Sea Ice Concentration Mask 
threshold = 0.80
meanSIarea = dsSI.SIarea.mean(dim='time')
Umask = grid.interp(meanSIarea,'X') > threshold
Vmask = grid.interp(meanSIarea,'Y') > threshold
Zmask = grid.interp(meanSIarea,['X','Y']) > threshold
Cmask = meanSIarea > threshold

# CHECK DATA FIRST

SIcontour = dsSI.SIarea.mean(dim='time')
levels =[1]
file = 'test'
Units = 'force'
ds = dsSI.isel(time=0)

sU, uEo, uEa, uE, uG, uRh, uC, uDer, uRHS, uLHS, uRES, uERR = decomp.Udecomp(ds,Units)
dsV,vEo, vEa, vE, vG, vRh, vC, vDer, vRHS, vLHS, vRES, vERR  = decomp.Vdecomp(ds,Units)
# decomp.mombudgetCheck(dsSI,uEo, uEa, uE, uG, uRh, uC, uDer, uRHS, uLHS, uRES, uERR,Units,file, SIcontour, levels)  
# decomp.mombudgetCheck(dsSI,vEo, vEa, vE, vG, vRh, vC, vDer, vRHS, vLHS, vRES, vERR,Units,file, SIcontour, levels)

