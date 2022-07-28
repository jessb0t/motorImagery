#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 09:08:17 2022

@author: kurtlehner
"""

from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import os
import pandas as pd
import seaborn as sns
import mne
from mne.baseline import rescale
from mne.stats import bootstrap_confidence_interval
from mne.io import concatenate_raws, read_raw_edf
from mne.time_frequency import tfr_multitaper
from mne.stats import permutation_cluster_1samp_test as pcluster_test
from mne_bids import BIDSPath, read_raw_bids
from mne.viz import plot_alignment, snapshot_brain_montage


# Load Matlab files 
ImageryData_pth = os.path.join(os.getcwd(),'sourcedata','imagery_basic', 'data')
ImageryGNRL_pth = os.path.join(os.getcwd(),'sourcedata','imagery_basic')

# Data load Function
def load_data(id,data_pth):
    imag = loadmat(os.path.join(data_pth,'{}_im_t_h.mat'.format(id)))
    overt = loadmat(os.path.join(data_pth,'{}_mot_t_h.mat'.format(id)))
    return imag,overt

def mne_conv(imag, overt, sampling_freq, scaling_f, n_channels):
    # Creating Info MNE Object
    mne_info= mne.create_info(n_channels, sampling_freq, ch_types=['ecog'] * n_channels)

    # Notch Filtering 
    filtered_data_Imag = []
    filtered_data_Overt = []
    for i in range(n_channels):
        channel_Imagdata = imag['data'][:,i] * scaling_f
        channel_Overtdata = overt['data'][:,i] * scaling_f

        filtered_Imag = mne.filter.notch_filter(channel_Imagdata, sampling_freq, np.arange(60,241,60))
        filtered_Overt = mne.filter.notch_filter(channel_Overtdata, sampling_freq, np.arange(60,241,60))

        filtered_data_Imag.append(filtered_Imag)
        filtered_data_Overt.append(filtered_Overt)

    # Creating MNE objects 
    mne_Overt = mne.io.RawArray(filtered_data_Overt, mne_info)
    mne_Imag = mne.io.RawArray(filtered_data_Imag, mne_info)
    
    # Stimulation events
    stim_events_overt = np.array([np.array([i, 0, stim[0]]) for i, stim in enumerate(overt['stim']) if overt['stim'][i-1]<stim]) # Stimulations
    stim_events_imag = np.array([np.array([i, 0, stim[0]]) for i, stim in enumerate(imag['stim']) if imag['stim'][i-1]<stim]) # Stimulations
    for ch in mne_Overt.ch_names:
        mne_Overt.add_events(stim_events_overt, stim_channel = ch)
        mne_Imag.add_events(stim_events_imag, stim_channel = ch)
    
    return mne_info, mne_Imag,mne_Overt

sampling_freq = 1000
scaling_f = 0.0298

#Load information
imag,overt=load_data('jc','sourcedata/imagery_basic/data')
n_channels = imag['data'].shape[1]
mne_Info, mne_Imag, mne_Overt=mne_conv(imag, overt, sampling_freq, scaling_f, n_channels)


#Set stim events
stim_events_overt = np.array([np.array([i, 0, stim[0]]) for i, stim in enumerate(overt['stim']) if overt['stim'][i-1]<stim]) # Stimulations
stim_events_imag = np.array([np.array([i, 0, stim[0]]) for i, stim in enumerate(imag['stim']) if imag['stim'][i-1]<stim]) # Stimulations

####

# let's explore some frequency bands
iter_freqs = [
    ('beta1', 12, 18),
    ('beta2', 18, 24),
    ('beta3', 24, 30),
    ('Gamma1', 30, 36),
    ('Gamma2', 36, 42),
    ('Gamma3', 70, 150)
]

# set epoching parameters
event_id, tmin, tmax = 12, -1., 3.
baseline = None


# get the header to extract events
raw_n = mne_Overt.copy()
stim_events = np.array([np.array([i, 0, stim[0]]) for i, stim in enumerate(overt['stim']) if overt['stim'][i-1]<stim]) # Stimulations


for ch in raw_n.ch_names:
    
    #Bandpass filter, epoch, envelope
    frequency_map = list()

    for band, fmin, fmax in iter_freqs:
        #(re)load the data to save memory
        raw_n = mne_Overt.copy()
        #raw_n.pick_types()  # we just look at gradiometers
        raw_n.load_data()
    
        # bandpass filter
        raw_n.filter(fmin, fmax, picks = [ch], n_jobs=1,  # use more jobs to speed up.
               l_trans_bandwidth=1,  # make sure filter params are the same
               h_trans_bandwidth=1)  # in each band and skip "auto" option.

        # epoch
        epochs = mne.Epochs(raw_n, events = stim_events, event_id = 11, tmin = tmin, tmax = tmax, baseline=None, picks = [ch],
                        preload=True)
        # remove evoked response
        # epochs.subtract_evoked()

        # get analytic signal (envelope)
        epochs.apply_hilbert(envelope=True, picks = [ch])
   
        frequency_map.append(((band, fmin, fmax), epochs.average()))
        del epochs
    del raw_n



#PSD plotting
#hand_epochs_Overt.plot_psd(fmin=12., fmax=30., average=True, picks='ecog', spatial_colors=False)
#hand_epochs_Img.plot_psd(fmin=12., fmax=30., average=True, picks='ecog', spatial_colors=False)

    # Helper function for plotting spread
    def stat_fun(x):
        """Return sum of squares."""
        return np.sum(x ** 2, axis=0)
    
    
    # Plot
    fig, axes = plt.subplots(6, 1, figsize=(10, 7), sharex=True, sharey=True)
    colors = plt.get_cmap('winter_r')(np.linspace(0, 1, 6))
       
    for ((freq_name, fmin, fmax), average), color, ax in zip(frequency_map, colors, axes.ravel()[::-1]):
        times = average.times * 1e3
        gfp = np.sum(average.data ** 2, axis=0)
        gfp = mne.baseline.rescale(gfp, times, baseline=(None, 0))
        ax.plot(times, gfp, label=freq_name, color=color, linewidth=2.5)
        ax.axhline(0, linestyle='--', color='grey', linewidth=2)
        ci_low, ci_up = bootstrap_confidence_interval(average.data, random_state=0,
                                                      stat_fun=stat_fun)
        ci_low = rescale(ci_low, average.times, baseline=(None, 0))
        ci_up = rescale(ci_up, average.times, baseline=(None, 0))
        ax.fill_between(times, gfp + ci_up, gfp - ci_low, color=color, alpha=0.3)
        ax.grid(True)
        ax.set_ylabel('GFP')
        ax.annotate('%s (%d-%dHz)' % (freq_name, fmin, fmax),
                            xy=(0.95, 0.8),
                            horizontalalignment='right',
                            xycoords='axes fraction')
        ax.set_xlim(-1000, 3000)
        ax.set_title(f'Channel {ch}')
    
    axes.ravel()[-1].set_xlabel('Time [ms]')




