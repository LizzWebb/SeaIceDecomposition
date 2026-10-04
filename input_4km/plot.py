# Plotting some of the data here

# import packages 
import xmitgcm
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path


# Import ds
data_dir = Path('/home/projects/def-straub/lizzwebb/mloch/arctic4km/run_mods')
ds       = xmitgcm.open_mdsdataset(data_dir,ref_date='1992-01-01 0:0:0', delta_t=240, swap_dims=False, \
                                   prefix = ['ds_hourly_snapshots'])
ds


