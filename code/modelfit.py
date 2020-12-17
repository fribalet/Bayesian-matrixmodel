#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 12:54:38 2020

@author: Paul Mattern
"""

import os
import pandas as pd
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt

_datafile = os.path.join(os.path.dirname(__file__), 'data', 'zinser_abundance_division_loss.nc')


def get_zinser_growth(correct_negative=False, minimum_sd=None, additive_sd=None):

    zinser_fig3 = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'Zinser_Figure3.csv'))

    tmp = zinser_fig3['Time of day'].values
    tmp -= tmp[0]
    time_hours = tmp.copy()
    delta = 0
    for i in range(1, len(tmp)):
        if tmp[i] < tmp[i-1]:
            delta += 24
        time_hours[i] += delta

    zinser_growth = zinser_fig3['PP'].values
    if correct_negative:
        zinser_growth[zinser_growth < 0] = 0.0

    sd = zinser_fig3['PP-sd'].values
    if additive_sd is not None:
        sd += additive_sd
    if minimum_sd is not None:
        sd = np.maximum(sd, minimum_sd)
    return zinser_growth, sd, time_hours


def get_zinser_abundance():
    with nc4.Dataset() as nc:
        return nc.variables['abundance'][:], nc.variables['abundance_sd'][:], nc.variables['abundance_time'][:]


def get_zinser_division(correct_negative=False):
    with nc4.Dataset(_datafile) as nc:
        division = nc.variables['division'][:]
        if correct_negative:
            division[division < 0] = 0.0
        return division, nc.variables['division_sd'][:], nc.variables['division_time'][:]


def get_zinser_carbonloss(correct_negative=False):
    with nc4.Dataset(_datafile) as nc:
        loss = nc.variables['carbonloss'][:]
        if correct_negative:
            loss[loss < 0] = 0.0
        return loss, nc.variables['carbonloss_sd'][:], nc.variables['carbonloss_time'][:]


def integrate(sample, time):
    dt = np.diff(time)
    return 0.5 * np.sum(dt * (sample[..., 1:] + sample[..., :-1]), axis=-1)


def extract_growth(nc_model, model_time_hours, compute_int=False):
    dt_h = model_time_hours[1] - model_time_hours[0]
    cell_count = nc_model.variables['cell_count'][:]
    if 'growth_vol_gain' in nc_model.variables:
        growth_size_gain = nc_model.variables['growth_vol_gain'][:]
    else:
        growth_size_gain = nc_model.variables['growth_size_gain'][:]
    # units: fg C cell^-1 h^-1
    growth = growth_size_gain/cell_count/dt_h
    if compute_int:
        return integrate(growth, model_time_hours), model_time_hours[-1] - model_time_hours[0]
    return growth, model_time_hours


def extract_loss(nc_model, model_time_hours, compute_int=False):
    dt_h = model_time_hours[1] - model_time_hours[0]
    cell_count = nc_model.variables['cell_count'][:]
    if 'resp_vol_loss' in nc_model.variables:
        resp_size_loss = nc_model.variables['resp_vol_loss'][:]
    else:
        resp_size_loss = nc_model.variables['resp_size_loss'][:]
    # units: fg C cell^-1 h^-1
    loss = resp_size_loss/cell_count/dt_h
    if compute_int:
        return integrate(loss, model_time_hours), model_time_hours[-1] - model_time_hours[0]
    return loss, model_time_hours


def compute_zinser_misfit(obstype, model_time_hours, growth_size_gain=None,
                          resp_size_loss=None, cell_count=None, nc_model=None,
                          misfit_type='chi2', return_plotdata=False,
                          create_debug_plot=False, debug_plot_title=None,
                          **extraargs):

    dt_h = model_time_hours[1] - model_time_hours[0]

    if obstype == 'growth':
        obs, obs_sd, obs_time_hours = get_zinser_growth(**extraargs)
    elif obstype == 'division':
        obs, obs_sd, obs_time_hours = get_zinser_division(**extraargs)
    elif obstype == 'carbonloss':
        obs, obs_sd, obs_time_hours = get_zinser_carbonloss(**extraargs)
    elif obstype == 'abundance':
        obs, obs_sd, obs_time_hours = get_zinser_abundance(**extraargs)
    else:
        raise ValueError('Input obstype must be one of "growth", "division", "carbonloss", or "abundance".')

    if cell_count is None:
        cell_count = nc_model.variables['cell_count'][:]

    if obstype == 'growth':
        if growth_size_gain is None:
            if 'growth_vol_gain' in nc_model.variables:
                growth_size_gain = nc_model.variables['growth_vol_gain'][:]
            else:
                growth_size_gain = nc_model.variables['growth_size_gain'][:]

        # units: fg C cell^-1 h^-1
        model = growth_size_gain/cell_count/dt_h
        model_time = model_time_hours
    elif obstype == 'division':

        model = (np.log(cell_count[:, 1:])-np.log(cell_count[:, :-1]))/(model_time_hours[1:]-model_time_hours[:-1])
        model_time = 0.5 * (model_time_hours[1:] + model_time_hours[:-1])
    elif obstype == 'carbonloss':
        if resp_size_loss is None:
            if 'resp_vol_loss' in nc_model.variables:
                resp_size_loss = nc_model.variables['resp_vol_loss'][:]
            else:
                resp_size_loss = nc_model.variables['resp_size_loss'][:]

        # units: fg C cell^-1 h^-1
        model = resp_size_loss/cell_count/dt_h
        model_time = model_time_hours
    elif obstype == 'abundance':
        model = cell_count * obs[0]
        model_time = model_time_hours

    model_interp = np.empty((model.shape[0], obs_time_hours.size))
    for isample in range(model.shape[0]):
        model_interp[isample, :] = np.interp(obs_time_hours, xp=model_time, fp=model[isample, :])

    if misfit_type == 'mean abs diff':
        misfit_elem = np.abs(model_interp - obs)
        misfit = np.mean(misfit_elem, axis=1)
    elif misfit_type == 'chi2':
        misfit_elem = ((model_interp - obs)/(obs_sd))**2
        misfit = np.sum(misfit_elem, axis=1)
    elif misfit_type == 'chi2/n':
        misfit_elem = ((model_interp - obs)/(obs_sd))**2
        misfit = np.mean(misfit_elem, axis=1)
    elif misfit_type == 'mean squared error':
        misfit_elem = (model_interp - obs)**2
        misfit = np.mean(misfit_elem, axis=1)
    else:
        raise ValueError('Unknown misfit_type "{}".'.format(misfit_type))

    if create_debug_plot:
        qq = np.percentile(model_interp, q=(5, 25, 50, 75, 95), axis=0)
        qq_misfit = np.percentile(misfit_elem, q=(5, 25, 50, 75, 95), axis=0)

        fig, axs = plt.subplots(nrows=2, sharex=True, figsize=(16, 9))

        ax = axs[0]
        ax.errorbar(obs_time_hours, obs, yerr=obs_sd, label='Zinser data', color='black')

        sc = ax.fill_between(obs_time_hours, qq[0, :], qq[-1, :], alpha=0.25)
        ax.fill_between(obs_time_hours, qq[1, :], qq[-2, :], alpha=0.5,
                        facecolor=sc.get_facecolor()[0])
        ax.plot(obs_time_hours, qq[2, :], color=sc.get_facecolor()[0][:-1],
                lw=2, label='model result')
        ax.grid(True)

        if debug_plot_title is not None:
            ax.set_title(debug_plot_title)

        ax = axs[1]
        sc = ax.fill_between(obs_time_hours, qq_misfit[0, :], qq_misfit[-1, :], alpha=0.25)
        ax.fill_between(obs_time_hours, qq_misfit[1, :], qq_misfit[-2, :], alpha=0.5, facecolor=sc.get_facecolor()[0])
        ax.plot(obs_time_hours, qq_misfit[2, :], color=sc.get_facecolor()[0][:-1], lw=2, label='model result')
        ax.grid(True)

    if return_plotdata:
        plotdata = {
            'model_quantiles': np.percentile(model, q=(5, 25, 50, 75, 95), axis=0),
            'model_time': model_time,
            'obs': obs,
            'obs_time': obs_time_hours,
            'obs_sd': obs_sd,
        }
        return misfit, plotdata
    return misfit


if __name__ == '__main__':

    datafile = 'local_files/data_exp_zs_20201105_g3_trackvol.nc'

    with nc4.Dataset(datafile) as nc:
        nc_model = nc['zinser/m4']
        model_time_hours = nc['zinser/time'][:] / 60
        compute_zinser_misfit('division',
                              model_time_hours=model_time_hours,
                              nc_model=nc_model)
