plot = 0
stream = 1

import time
import mne
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import filtfilt, butter, windows
from pylsl import StreamInfo, StreamOutlet, local_clock

# Bandpass filter
def bandpass_filter(signal, low_freq, high_freq, sfreq, order=5):
    nyq = 0.5 * sfreq
    low = low_freq / nyq
    high = high_freq / nyq
    b, a = butter(order, [low, high], btype='band')
    return filtfilt(b, a, signal)

# Generate pink noise
def generate_pink_noise(duration_s, sfreq, num_channels, random_seed=6):
    #np.random.seed(random_seed)
    white_noise = np.random.randn(num_channels, int(duration_s * sfreq))
    pink_noise = np.apply_along_axis(lambda e: np.convolve(e, np.ones((100,))/100, mode='same'), axis=1, arr=white_noise)
    # Approximation of pink noise:
    # A simple moving average filter (SMA; convolution of a filter kernel (=sequence of weights) with the white_noise)
    # is similar to a low-pass filter, effectively reducing high-frequency components, which is why low-frequency components in the
    # PSD have higher power than high-frequency components, resulting in a 1/f characteristic of pink noise, similar to EEG data. 
    return pink_noise

# Simulate EEG data with beta-ERD:
# 1) Generate pink noise
# 2) Bandpass filter beta band (13-20 Hz)
# 3) Scale amplitude of filtered beta between 1-4s to simulate ERD
# 4) Remove beta from pink noise signal
# 5) Add scaled beta
# Trial: 1s baseline, 3s beta-ERD, 5s inter-trial interval
def simulate_trial_with_beta_erd(sfreq, n_channels, ch_names, duration_s, erd_start_s, erd_duration_s):
    # Generate pink noise as base EEG signal
    eeg_data = generate_pink_noise(duration_s, sfreq, n_channels)

    # Index for C3 channel
    c3_index = ch_names.index("C3")

    # Isolate the beta band component from the C3 channel
    beta_signal = bandpass_filter(eeg_data[c3_index, :], 13, 20, sfreq)

    # Scaling factor for beta-ERD (reduce beta amplitude to simulate ERD)
    reduction_factor = np.sqrt(0.1)  # Reduce power by ...

    # Scale down the beta signal for the duration of the ERD
    start_sample = int(erd_start_s * sfreq)
    end_sample = int((erd_start_s + erd_duration_s) * sfreq)
    beta_signal[start_sample:end_sample] *= reduction_factor

    # Remove the original beta component 
    eeg_data[c3_index, :] -= bandpass_filter(eeg_data[c3_index, :], 13, 20, sfreq)
        
    # Add the scaled beta component
    eeg_data[c3_index, :] += beta_signal
    eeg_data[c3_index, :] = beta_signal

    # Increase the magnitude of the EEG data
    eeg_data = eeg_data * 100
    
    return eeg_data

def simulate_eeg_data_with_variable_erd(sfreq, n_channels, ch_names, erd_durations, n_trials=100):
    total_duration_s = n_trials * 10  # 10 seconds per trial for 100 trials
    all_trials_data = np.empty((n_channels, 0))

    for trial in range(n_trials):
        # Use the ERD duration for the current trial
        erd_duration_s = erd_durations[trial]
        # iti_duration_s = 10 - 1 - erd_duration_s

        # Generate trial data with the selected ERD duration
        trial_data = simulate_trial_with_beta_erd(sfreq, n_channels, ch_names, 
                                                     duration_s=10, erd_start_s=1, erd_duration_s=erd_duration_s)
        all_trials_data = np.concatenate((all_trials_data, trial_data), axis=1)
    
    return all_trials_data

# Settings
sfreq = 1000  # Sampling rate
n_channels = 2  # Number of EEG channels
ch_names = ['C3', 'C16']
# n_channels = 16  # Number of EEG channels
# ch_names = ['C3', 'C16', 'C29', 'D7', 'D4', 'C4', 'C7', 'D23', 'B26', 'D19', 'B22', 'D31', 'A8', 'B5', 'B11', 'A15']

erd_durations = np.round(np.random.uniform(2, 5, size=100), 2)

# Simulate EEG data with beta-ERD
simulated_data = simulate_eeg_data_with_variable_erd(sfreq, n_channels, ch_names, erd_durations)

# Create MNE Raw object for simulated data
info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types=['eeg'] * n_channels)
raw_simulated = mne.io.RawArray(simulated_data, info)

if stream:
    # Create EEG stream info and outlet
    eeg_info = StreamInfo('sim_ERD', 'EEG', n_channels, sfreq, 'float32', 'myuid_eeg')
    chns = eeg_info.desc().append_child("channels")
    for chan_ix, label in enumerate(ch_names):
        ch = chns.append_child("channel")
        ch.append_child_value("label", label)
        ch.append_child_value("unit", "microvolts")
        ch.append_child_value("type", "EEG")
    eeg_outlet = StreamOutlet(eeg_info)

    # Create Marker stream info and outlet
    mar_info = StreamInfo('sim_Markers', 'Markers', 1, 0, 'string', 'myuid_beh')
    mar_outlet = StreamOutlet(mar_info)

    print("Now sending data...")
    start_time = local_clock()
    sent_samples = 0
    next_trial_start_time = start_time  # Initialize next trial start time

    for trial in range(len(erd_durations)):
        trial_start_time = next_trial_start_time  # Mark the start of the current trial
        trial_start_timestamp = f"{trial_start_time - start_time:.2f}"
        # baseline_marker = f"bsl:{trial_start_time - start_time:.2f}"
        # mar_outlet.push_sample([baseline_marker])
        # print(baseline_marker)

        # ERD Marker should be sent 1 second after the baseline marker
        erd_onset_time = trial_start_time + 1
        erd_onset_timestamp = f"{erd_onset_time - start_time:.2f}"
        # ITI Marker should be sent immediately after the ERD duration ends
        iti_onset_time = erd_onset_time + erd_durations[trial]
        iti_onset_timestamp = f"{iti_onset_time - start_time:.2f}"
        # Next trial's baseline marker should be exactly 10 seconds after the current trial's baseline marker
        next_trial_start_time = trial_start_time + 10

        while True:
            current_time = local_clock()
            if current_time >= erd_onset_time and current_time < iti_onset_time:
                erd_marker = f"erd"
                mar_outlet.push_sample([erd_marker])
                print(erd_marker)
                erd_onset_time = float('inf')  # Prevent repeated sending

            if current_time >= iti_onset_time:
                iti_marker = f"iti"
                mar_outlet.push_sample([iti_marker])
                print(iti_marker)
                iti_onset_time = float('inf')  # Prevent repeated sending

            # Break out of the loop once it's time to start the next trial
            if current_time >= next_trial_start_time:
                break

            # Send EEG data
            required_samples = int(sfreq * (current_time - start_time)) - sent_samples
            if required_samples > 0:
                data = raw_simulated.get_data(start=sent_samples, stop=sent_samples + required_samples, return_times=False)
                eeg_outlet.push_chunk(data.T.tolist())
                sent_samples += required_samples

            # Sleep for a short duration to control the rate of data sending
            time.sleep(1/sfreq)

if plot:
    # Plot the simulated EEG signal
    raw_simulated.plot(n_channels=n_channels, scalings='auto', title='Simulated EEG Signal')
    plt.show()

    # Extract C3 data from simulated EEG
    c3_data = raw_simulated.get_data(picks='C3')

    # Define the start and stop sample indices for each period
    baseline_start = 0
    baseline_stop = int(sfreq * 1)  # 1 second of baseline
    beta_erd_start = baseline_stop
    beta_erd_stop = int(beta_erd_start + sfreq * 3)  # 3 seconds of beta-ERD
    iti_start = beta_erd_stop
    iti_stop = int(iti_start + sfreq * 5)  # 5 seconds of inter-trial interval

    # Plot PSD for the three periods
    fig, axs = plt.subplots(3, 1, figsize=(10, 7), sharex=True, sharey=True)

    # Plot PSD for the three periods with consistent windowing and overlap
    def plot_psd(data, start, stop, ax, title):
        # Extract the segment while keeping the data 2D
        segment_data = data[:, start:stop]
        segment_info = mne.create_info(ch_names=['C3'], sfreq=sfreq, ch_types=['eeg'])
        segment_raw = mne.io.RawArray(segment_data, segment_info)
        
        # Define window and overlap
        n_per_seg = int(sfreq / 2)  # half-second segments
        window = windows.hann(n_per_seg)
        overlap = n_per_seg // 2
        
        # Plot the power spectral density
        segment_raw.plot_psd(fmin=1, fmax=40, ax=ax, spatial_colors=False, 
                            n_per_seg=n_per_seg, n_overlap=overlap, window=window, show=False)
        ax.set_title(title)

    # Plot the PSD for each period with consistent windowing and overlap
    plot_psd(c3_data, baseline_start, baseline_stop, axs[0], 'Baseline')
    plot_psd(c3_data, beta_erd_start, beta_erd_stop, axs[1], 'Beta-ERD')
    plot_psd(c3_data, iti_start, iti_stop, axs[2], 'Inter-Trial Interval')

    plt.tight_layout()
    plt.show()



