from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import os
import numpy as np
import sys 
import mne 
# from nilearn import plotting
# from nimare import utils

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

def filter_frequency_bands(mne_Imag, mne_Overt, freq_limit, picks=[]):
    if len(picks)==0:
        picks = mne_Overt.ch_names
    freq_mne_Overt = mne_Overt.copy().filter(freq_limit[0], freq_limit[1], picks=picks).pick(picks)
    freq_mne_Imag = mne_Imag.copy().filter(freq_limit[0], freq_limit[1], picks=picks).pick(picks)
    return freq_mne_Overt, freq_mne_Imag

sampling_freq = 1000
scaling_f = 0.0298

#beta1
beta1 = (12, 18)
#beta2
beta2 = (18, 24)
#beta3
beta3 = (24, 30)
#gamma1
gamma1 =(30, 36)
#gamma2
gamma2 =(36, 42)
#gammaH
gammaH =(70, 150)

imag,overt=load_data('jc','sourcedata/imagery_basic/data')
n_channels = imag['data'].shape[1]
mne_Imag,mne_Overt=mne_conv(imag, overt, sampling_freq, scaling_f, n_channels)
stim_events = np.array([np.array([i, 0, stim[0]]) for i, stim in enumerate(overt['stim']) if overt['stim'][i-1]<stim]) # Stimulations
for mne_obj in [mne_Imag,mne_Overt]:
    for ch in mne_obj.ch_names:
        mne_obj.add_events(stim_events, stim_channel = ch)

picks = [19, 20, 21, 26, 27, 28, 35, 36]
freq_mne_Overt, freq_mne_Imag = filter_frequency_bands(mne_Imag, mne_Overt, beta2)
# freq_mne_Overt.plot(events=stim_events, event_color={11:'g', 12:'b', -1:'w'})

# apply hilbert transform
hilbert_freq_mne_Overt = freq_mne_Overt.copy().apply_hilbert(picks = picks, envelope=True)#.pick(picks)
# hilbert_freq_mne_Overt.plot(events=stim_events, event_color={11:'g', 12:'b', -1:'w'})

# apply low sustained and dynamic filter
sustained_hilbert_freq_mne_Overt = hilbert_freq_mne_Overt.copy().filter(None, 0.4, picks=picks)
sustained_hilbert_freq_mne_Overt.pick(picks).plot(events=stim_events, event_color={11:'g', 12:'b', -1:'w'})

dynamic_hilbert_freq_mne_Overt = hilbert_freq_mne_Overt.copy().filter(0.4,None, picks=picks)
dynamic_hilbert_freq_mne_Overt.plot(events=stim_events, event_color={11:'g', 12:'b', -1:'w'})





