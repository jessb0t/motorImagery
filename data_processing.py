from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import os
import numpy as np
import sys 
import mne
# ------------------------Functions--------------------------------------
# Data load Function
def load_data(id,data_pth):
    imag = loadmat(os.path.join(data_pth,'{}_im_t_h.mat'.format(id)))
    overt = loadmat(os.path.join(data_pth,'{}_mot_t_h.mat'.format(id)))
    return imag,overt

def mne_conv(imag, overt, sampling_freq, scaling_f, n_channels):
    # Creating Info MNE Object
    mne_info= mne.create_info(n_channels, sampling_freq)#, ch_types=['ecog']*n_channels)

    # Notch Filtering 
    filtered_data_Imag = []
    filtered_data_Overt = []
    for i in range(n_channels):
        channel_Imagdata = imag['data'][:,i] * scaling_f
        channel_Overtdata = overt['data'][520:,i] * scaling_f

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

def filter_frequency_bands(mne_Imag, mne_Overt, freq_limit, picks=[]):
    if len(picks)==0:
        picks = mne_Overt.ch_names
    freq_mne_Overt = mne_Overt.copy().filter(freq_limit[0], freq_limit[1], picks=picks).pick(picks)
    freq_mne_Imag = mne_Imag.copy().filter(freq_limit[0], freq_limit[1], picks=picks).pick(picks)
    return freq_mne_Overt, freq_mne_Imag

def get_epochs_data(hilbert_freq_mne, event_id):
    epochs_mne_data = mne.Epochs(hilbert_freq_mne.copy(), events=stim_events, tmin=-1, tmax=3)[event_id]
    return epochs_mne_data

# -----------------variables--------------------
freq_bands = { 
    'beta1' : (12, 18),
    'beta2' : (18, 24),
    'beta3' : (24, 30),
    'gamma1' :(30, 36),
    'gamma2' :(36, 42),
    'gammaH' :(70, 150)
}
sampling_freq = 1000
scaling_f = 0.0298
data_pth = os.path.join(os.getcwd(),'sourcedata','imagery_basic','data')
imag,overt=load_data('jc',data_pth)
n_channels = imag['data'].shape[1]
mne_Info, mne_Imag ,mne_Overt=mne_conv(imag, overt, sampling_freq, scaling_f, n_channels)

stim_events = np.array([np.array([i, 0, stim[0]]) for i, stim in enumerate(overt['stim']) if overt['stim'][i-1]<stim]) # Stimulations
for mne_obj in [mne_Imag,mne_Overt]:
    for ch in mne_obj.ch_names:
        mne_obj.add_events(stim_events, stim_channel = ch)

# overt['stim'][-520:].max()

save_path = os.path.join(os.getcwd(),'processed_data')

hilbert_freq_mne_Overt = dict()
hilbert_freq_mne_Imag = dict()
hand_epochs_Img = dict()
hand_epochs_Overt = dict()
tng_epochs_Img = dict()
tng_epochs_Overt = dict()

hand_epochs_overt_data = dict()
hand_epochs_img_data = dict()
hilbert_freq_img_data = dict()
hilbert_freq_overt_data = dict()
picks = [19, 20, 21, 26, 27, 28, 35, 36]

for band in freq_bands:
    freq_mne_Overt, freq_mne_Imag = filter_frequency_bands(mne_Imag, mne_Overt, freq_bands[band])
    hilbert_freq_mne_Overt[band] = freq_mne_Overt.copy().apply_hilbert(picks = picks, envelope=True)
    hilbert_freq_mne_Imag[band] = freq_mne_Imag.copy().apply_hilbert(picks = picks, envelope=True)
    
    hand_epochs_Img[band] = get_epochs_data(hilbert_freq_mne_Imag[band].copy(), '12')
    hand_epochs_Overt[band] = get_epochs_data(hilbert_freq_mne_Overt[band].copy(), '12')
    
    tng_epochs_Img[band] = get_epochs_data(hilbert_freq_mne_Imag[band].copy(), '11')
    tng_epochs_Overt[band] = get_epochs_data(hilbert_freq_mne_Overt[band].copy(), '11')
    
    save_path = os.path.join(os.getcwd(),'processed_data')

    # np.save(os.path.join(save_path,'hand_Epochs_Img_{}.npy'.format(band)),np.array(hand_epochs_Img[band]._get_data()))
    # np.save(os.path.join(save_path,'hand_Epochs_Overt_{}.npy'.format(band)),np.array(hand_epochs_Overt[band]._get_data()))
    # np.save(os.path.join(save_path,'Tng_Epochs_Overt_{}.npy'.format(band)),np.array(tng_epochs_Overt[band]._get_data()))
    # np.save(os.path.join(save_path,'Tng_Epochs_Img_{}.npy'.format(band)),np.array(tng_epochs_Img[band]._get_data()))
    
    hand_epochs_overt_data[band] = hand_epochs_Overt[band]._get_data()
    hand_epochs_img_data[band] = hand_epochs_Img[band]._get_data()
    
    hilbert_freq_img_data[band] = hilbert_freq_mne_Imag[band].get_data()
    hilbert_freq_overt_data[band] = hilbert_freq_mne_Overt[band].get_data()

    print('saved '+ band)










