#!/usr/bin/env python3

import pandas as pd
import numpy as np
import netCDF4 as nc4
import sys

zinser_data = pd.read_csv('Zinser_Figure2A.csv')

#zinser_counts = zinser.values[:,1].astype(int) # cells A column
#zinser_counts = zinser.values[:,2].astype(int) # cells B column
zinser_abundances = np.mean(zinser_data.values, axis=1).astype(int) # mean of both columns

if len(sys.argv) > 1:
    with nc4.Dataset(sys.argv[1], 'a') as nc:
        # check time
        assert nc.variables['time'].size == len(zinser_abundances)
        nc.createVariable('abundance', int, ('time',), fill_value=False)
        nc.variables['abundance'][:] = zinser_abundances
        nc.variables['abundance'].units = 'cells ml-1'
        nc.variables['abundance'].long_name = 'cell abundance'
