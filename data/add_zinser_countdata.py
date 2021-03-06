#!/usr/bin/env python3

import pandas as pd
import numpy as np
import netCDF4 as nc4
import sys

zinser_data = pd.read_csv('Zinser_Figure2A.csv')
zinser_abundances = (0.5 * (zinser_data['cells A'].values + zinser_data['cells B'].values)).astype(int) # mean of both columns

if len(sys.argv) > 1:
    with nc4.Dataset(sys.argv[1], 'a') as nc:
        # check time
        assert nc.variables['time'].size == len(zinser_abundances)

        if 'abundance' not in nc.variables:
            nc.createVariable('abundance', int, ('time',), fill_value=False)
        nc.variables['abundance'][:] = zinser_abundances
        nc.variables['abundance'].units = 'cells ml-1'
        nc.variables['abundance'].long_name = 'cell abundance'
        
        # For Zinser data, 20 microliter of sample (diluted 1/100) was analyzed, so divide the cell abundance by the volume to determine how many cells were counted by the instrument.
        if 'count' not in nc.variables:
            nc.createVariable('count', int, ('time',), fill_value=False)
        nc.variables['count'][:] = (np.round(zinser_abundances * 0.0002)).astype(int) # 20 microliter = 0.02 milliliter 
        nc.variables['count'].units = 'cells'
        nc.variables['count'].long_name = 'cell count'
