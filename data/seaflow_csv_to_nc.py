#!/usr/bin/env python3

import netCDF4 as nc4
import numpy as np
import os
from dateutil.parser import parse
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

description = '''This script regrids the data in a SeaFlow csv-file into a format
compatible with the matrix population model. The parameters v_min,
m and delta_v_inv are set below, m and delta_v_inv can also be 
specified as input parameters. If --createplots is active, both
the old and new size grid are plotted along with the old and regridded 
data set. The regridded data is saved in NetCDF format.'''

def regrid(v0, w0, v1, permit_smallgrid=False):
    # assert that old grid is contained in new one
    if  v1[0] > v0[0]:
        raise ValueError('Value of v_min ({}) is too large to capture data (must be <= {}).'.format(v1[0], v0[0]))
    if not permit_smallgrid:
        assert v1[-1] >= v0[-1]
    # check sizes
    if v0.size != w0.shape[0]+1:
        raise ValueError('Size of v0 ({}) does not match shape of w0 {}.'.format(v0.size, w0.shape))

    w1 = np.zeros((v1.size-1,w0.shape[1]))
    #v_data = data['v_edges']

    i1 = 1 # first right edge
    for i0 in range(1,v0.size):
        if v0[i0] < v1[i1]:
            if v0[i0-1] >= v1[i1-1]:
                w1[i1-1,:] += w0[i0-1,:]
            else:
                raise RuntimeError('This should not happen.')
        elif i1+1 == v1.size:
            # can only happen for permit_smallgrid
            a = (v1[i1]-v0[i0-1])/(v0[i0]-v0[i0-1])
            # left side
            w1[i1-1,:] += w0[i0-1,:]*a
            # right side is not covered by v1
        elif v0[i0] < v1[i1+1]:
            a = (v1[i1]-v0[i0-1])/(v0[i0]-v0[i0-1])
            # left side
            w1[i1-1,:] += w0[i0-1,:]*a
            # right side
            w1[i1,:] += w0[i0-1,:]*(1-a)
            i1 += 1
        else:
            raise NotImplementedError('Case not yet covered, chose a v_min smaller or equal to smallest value in data grid.')
    
    if v1[-1] >= v0[-1]:
        assert np.all(np.abs(np.sum(w0,axis=0)-np.sum(w1,axis=0))<1e-9)
    else:
        print('maximum loss of data due to reduction in grid coverage: {}'.format(np.max(np.abs(np.sum(w0,axis=0)-np.sum(w1,axis=0)))))
    return w1


def add_colorbar(ax, **cbarargs):
    axins_cbar = inset_axes(ax, width='3%', height='90%', loc=5, bbox_to_anchor=(0.05,0.0,1,1), bbox_transform=ax.transAxes)
    mpl.colorbar.ColorbarBase(axins_cbar, orientation='vertical', **cbarargs)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('csvfile', type=str, help='The input csv-file.')
    parser.add_argument('--vmin', type=float, default=0.0135, help='The value of the parameter v_min.')
    parser.add_argument('--m', type=int, default=25, help='The value of the parameter m.')
    parser.add_argument('--deltavinv', type=int, default=8, help='The value of the parameter delta_v_inv.')
    parser.add_argument('--output', type=str, default=None, help='The name of the output file.')
    parser.add_argument('--minbinindex', type=int, default=5, help='A zero-based index denoting the start of bins in the data.')
    parser.add_argument('--createplots', action='store_true', help='Show plots of the data.')
    parser.add_argument('--nosmallgrid', action='store_true', help='Do not permit the upper bounds of the data not to be contained by the model grid.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing file if output file exists.')
    parser.add_argument('--startday', type=int, default=None, help='The first day to include in output file (default: include every day).')
    parser.add_argument('--endday', type=int, default=None, help='The last day to start include in output file (default: include every day).')
    
    args = parser.parse_args()

    data_raw = pd.read_csv(args.csvfile)
    
    data = {}
    data['v_centers'] = np.array([float(x) for x in data_raw.columns[args.minbinindex:]])
    data['v_edges'] = np.empty(data['v_centers'].size + 1)
    data['v_edges'][1:-1] = 0.5 * (data['v_centers'][1:] + data['v_centers'][:-1])
    data['v_edges'][0] = 0.0001*np.round(10000 * (data['v_centers'][0] - (data['v_edges'][1] - data['v_centers'][0])))
    data['v_edges'][-1] = data['v_centers'][-1] + (data['v_centers'][-1] - data['v_edges'][-2])

    date = [parse(d) for d in data_raw['date'].values]
    data['t_min'] = np.array([(d-date[0]).total_seconds()//60 for d in date])
    
    filesuffix = ''
    if args.startday is not None or args.endday is not None:
        day = data['t_min']//1400 + 1
        index = np.ones_like(day, dtype=bool)

        if args.startday is not None:
            index &= day >= args.startday
        if args.endday is not None:
            index &= day <= args.endday
        
        data_raw = data_raw.iloc[index,:]
        data['t_min'] = data['t_min'][index]
        data['t_min'] -= data['t_min'][0]
        

        if args.startday is None:
            args.startday = day[0]
        if args.endday is None:
            args.endday = day[-1]
        if args.startday == args.endday:
            filesuffix = '_day{}'.format(args.startday)
        else:
            filesuffix = '_day{}to{}'.format(args.startday, args.endday)

    data['data'] = data_raw.iloc[:,args.minbinindex:].values
    data['par'] = data_raw['PAR'].values


    #
    # perform regridding
    #

    v_min = args.vmin
    m = args.m
    delta_v_inv = args.deltavinv

    if args.output is None:
        args.output = args.csvfile.replace('.csv','') + '_regrid-{}-{}{}.nc'.format(m, delta_v_inv, filesuffix)
    if os.path.isfile(args.output) and not args.overwrite:
        raise RuntimeError('File "{}" already exists (use --overwrite to overwrite it).'.format(args.output))
    
    delta_v = 1.0/delta_v_inv
    v = v_min * 2**(np.arange(m+1)*delta_v) # to get m intervals, we need m+1 edges

    if args.createplots:
        fig,ax = plt.subplots()
        ax.plot(data['v_edges'], np.ones_like(data['v_edges']), marker='s', label='data v edges')
        ax.plot(v, np.zeros_like(v), marker='o', label='model v edges')
        ax.legend()
        ax.set(ylim=(-1,2), yticks=[])
        #plt.show()

    # data['data'] is num_t x num_v
    w = regrid(data['v_edges'], data['data'].T, v, permit_smallgrid=not args.nosmallgrid)

    if args.createplots:
        fig,axs = plt.subplots(3,1,sharex=True,figsize=(12,10))
        ax = axs[0]
        ax.plot(data['t_min'], data['par'], color='gold')
        ax.set(ylabel='PAR (µmol photons m$^{-2}$ s$^{-1}$)')
        ax = axs[1]
        pc = ax.pcolormesh(data['t_min'], data['v_edges'], data['data'].T/np.sum(data['data'],axis=1), rasterized=True)
        add_colorbar(ax, label='size class proportion', norm=pc.norm, cmap=pc.cmap)
        ax.set(ylabel='size ($\mu$m$^3$)')
        ax = axs[2]
        pc = ax.pcolormesh(data['t_min'], v, w/np.sum(w,axis=0), rasterized=True)
        add_colorbar(ax, label='size class proportion', norm=pc.norm, cmap=pc.cmap)
        ax.set(xlabel='minutes since start', ylabel='size ($\mu$m$^3$)', ylim=axs[1].get_ylim())
        ax.text(0.01, 0.98, 'm = {}, delta_v_inv = {}'.format(m, delta_v_inv), ha='left', va='top', transform=ax.transAxes, size=16)
        
        #fig.patch.set_alpha(0.0)
        #fname = 'SeaFlow_SizeDist_regrid-{}-{}.pdf'.format(m, delta_v_inv)
        #fig.savefig(fname, bbox_inches='tight')

    #
    # write to netCDF
    #
    
    print('writing to "{}"'.format(args.output))
    with nc4.Dataset(args.output, 'w') as nc:
        nc.createDimension('time', w.shape[1])
        nc.createDimension('size', w.shape[0])
        nc.createDimension('size_bounds', w.shape[0]+1)

        nc.createVariable('time', int, ('time',), fill_value=False)
        nc.variables['time'][:] = data['t_min']
        nc.variables['time'].units = 'minutes since start of experiment'

        nc.createVariable('size', float, ('size',), fill_value=False)
        nc.variables['size'][:] = v[:-1]
        nc.variables['size'].units = 'um3'
        nc.variables['size'].long_name = 'the lower bound for each size class (used as the size of each size class in the model)'

        nc.createVariable('size_bounds', float, ('size_bounds',), fill_value=False)
        nc.variables['size_bounds'][:] = v
        nc.variables['size_bounds'].units = 'um3'
        nc.variables['size_bounds'].long_name = 'the lower and upper bounds for each size class'

        nc.createVariable('w_obs', float, ('size','time'), fill_value=False)
        nc.variables['w_obs'][:] = w/np.sum(w,axis=0)
        nc.variables['w_obs'].units = 'proportion'
        nc.variables['w_obs'].long_name = 'proportion of the population in each size class'

        nc.createVariable('PAR', float, ('time',), fill_value=False)
        nc.variables['PAR'][:] = data['par']
        nc.variables['PAR'].units = 'umol photons/m2/s'

        nc.createVariable('m', int, fill_value=False)
        nc.variables['m'][:] = m
        nc.variables['m'].long_name = 'the number of size classes'

        nc.createVariable('delta_v_inv', int, fill_value=False)
        nc.variables['delta_v_inv'][:] = delta_v_inv
        nc.variables['delta_v_inv'].long_name = 'the inverse of delta_v, spacifying the logarithmic scaling of size classes'

        nc.createVariable('v_min', float, fill_value=False)
        nc.variables['v_min'][:] = v_min
        nc.variables['v_min'].long_name = 'the lower bound of the smallest size class (used as the size of the smallest size class in the model)'

    if args.createplots:
        plt.show()
