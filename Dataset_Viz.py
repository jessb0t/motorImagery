from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import os 
import sys 
import mne 
from nilearn import plotting
from nimare import utils

# Load Matlab files 
ImageryData_pth = os.path.join(os.getcwd(),'sourcedata','imagery_basic', 'data')
ImageryGNRL_pth = os.path.join(os.getcwd(),'sourcedata','imagery_basic')

# Data load Function
def load_data(id,data_pth):
    imag = loadmat(os.path.join(data_pth,'{}_im_t_h.mat'.format(id)))
    overt = loadmat(os.path.join(data_pth,'{}_mot_t_h.mat'.format(id)))
    return imag,overt


# # Participant "bp"
# bp_imag = loadmat(os.path.join(ImageryData_pth,'bp_im_t_h.mat'))
# bp_overt = loadmat(os.path.join(ImageryData_pth,'bp_mot_t_h.mat'))

# # Participant "fp"
# fp_imag = loadmat(os.path.join(ImageryData_pth,'fp_im_t_h.mat'))
# fp_overt = loadmat(os.path.join(ImageryData_pth,'fp_mot_t_h.mat'))

# # Participant "hh"
# hh_imag = loadmat(os.path.join(ImageryData_pth,'hh_im_t_h.mat'))
# hh_overt = loadmat(os.path.join(ImageryData_pth,'hh_mot_t_h.mat'))

# # Participant "jc"
# jc_imag = loadmat(os.path.join(ImageryData_pth,'jc_im_t_h.mat'))
# jc_overt = loadmat(os.path.join(ImageryData_pth,'jc_mot_t_h.mat'))

# # Participant "jm"
# jm_imag = loadmat(os.path.join(ImageryData_pth,'jm_im_t_h.mat'))
# jm_overt = loadmat(os.path.join(ImageryData_pth,'jm_mot_t_h.mat'))

# # Participant "rh"
# rh_imag = loadmat(os.path.join(ImageryData_pth,'rh_im_t_h.mat'))
# rh_overt = loadmat(os.path.join(ImageryData_pth,'rh_mot_t_h.mat'))

# # Participant "rr"
# rr_imag = loadmat(os.path.join(ImageryData_pth,'rr_im_t_h.mat'))
# rr_overt = loadmat(os.path.join(ImageryData_pth,'rr_mot_t_h.mat'))

def mne_conv(imag,overt):
    # Creating Info MNE Object
    n_channels = imag.shape[1]
    sampling_freq = 1000

    scaling_f = 0.0298

    mne_info= mne.create_info(n_channels, sampling_freq)

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
    
    return mne_Imag,mne_Overt

 




    

# Plotting 
# Data
hand_events = [(i, 0, stim[0]) for i, stim in enumerate(bp_overt['stim']) if bp_overt['stim'][i-1]<stim and stim==12] # Hand Stimulations
stim_events = np.array([np.array([i, 0, stim[0]]) for i, stim in enumerate(bp_overt['stim']) if bp_overt['stim'][i-1]<stim]) # Stimulations
jc_stim_events = np.array([np.array([i, 0, stim[0]]) for i, stim in enumerate(bp_overt['stim']) if bp_overt['stim'][i-1]<stim]) # Stimulations
for ch in mne_bpOvert.ch_names:
    mne_bpOvert.add_events(stim_events, stim_channel = ch)

mne_bpOvert.plot(n_channels = n_channels, events = stim_events, event_color={11:'g', 12:'b', -1:'w'})

# Electrode Location
# load loc file
loc_pth = os.path.join(ImageryGNRL_pth, 'locs','bp_electrodes.mat')
bp_locs = loadmat(loc_pth) 
plt.figure(figsize=(8, 8))
locs = bp_locs['electrodes']
view = plotting.view_markers(utils.tal2mni(locs),
                             marker_labels=['%d'%k for k in np.arange(locs.shape[0])],
                             marker_color='purple',
                             marker_size=5)
view.open_in_browser()

# Pre-processing
# Filtering relevant frequency bands
picks_bp = [10, 11, 18, 19, 20, 25, 26, 27]
picks_jc = [19, 20, 21, 26, 27, 28, 35, 36]
#beta1
bp_beta1 = mne_bpOvert.copy().filter(12, 18, picks=picks).pick(picks)

#beta2
bp_beta2 = mne_bpOvert.copy().filter(18, 24, picks=picks).pick(picks)

#beta3
bp_beta3 = mne_bpOvert.copy().filter(24, 30, picks=picks).pick(picks)

#gamma1
bp_gamma1 = mne_bpOvert.copy().filter(30, 36, picks=picks).pick(picks)

#gamma2
bp_gamma2 = mne_bpOvert.copy().filter(36, 42, picks=picks).pick(picks)

#gammaH
bp_gammaH = mne_bpOvert.copy().filter(70, 150, picks=picks_jc).pick(picks_jc)

bp_beta1.plot(events = stim_events, event_color={11:'g', 12:'b', -1:'w'})

print('plotting')










