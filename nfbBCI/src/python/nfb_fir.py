from pylsl import StreamInlet, StreamInfo, StreamOutlet, resolve_stream
import numpy as np
from scipy.signal import firwin, lfilter, group_delay
import requests
import time

# Replace with Arduino's IP address
vibFeedback = 0
bsl = 10  # s
url = 'http://192.168.137.223/'

def toggle_vib(state):
    if state.lower() == 'on':
        requests.get(url + 'vibMotor=ON')
    elif state.lower() == 'off':
        requests.get(url + 'vibMotor=OFF')
    else:
        print("Invalid state. Use 'on' or 'off'.")

class Buffer:
    def __init__(self, size, num_channels):
        self.size = size
        self.num_channels = int(num_channels)
        self.data = np.zeros((num_channels, size), dtype=float)

    def add(self, samples):
        samples = np.atleast_2d(samples)
        num_new_samples = samples.shape[1]
        if num_new_samples >= self.size:
            self.data = samples[:, -self.size:]
        else:
            self.data = np.roll(self.data, -num_new_samples, axis=1)
            self.data[:, -num_new_samples:] = samples

    def get_latest_window(self):
        return self.data[:, -self.size:]
    
    def get_data(self):
        return self.data

# FIR Bandpass filter with group delay calculation
def bandpass_filter(signal, low_freq, high_freq, sfreq, numtaps=101):
    # Design FIR filter
    nyq = 0.5 * sfreq
    low = low_freq / nyq
    high = high_freq / nyq
    taps = firwin(numtaps, [low, high], pass_zero=False)
    
    # Compute group delay
    w, gd = group_delay((taps, 1))
    
    # Calculate the average group delay in samples
    avg_gd = int(np.mean(gd))
    # print(avg_gd)
    
    # Apply FIR filter in one direction
    filtered_signal = lfilter(taps, 1.0, signal, axis=1)
    
    return filtered_signal, avg_gd

def calculate_power(data, instantaneous=False):
    if instantaneous:
        # Instantaneous power calculation
        power = data ** 2
        return power
    else:
        # Median power calculation
        median_power = np.median(data ** 2, axis=1)
        return median_power

# LSL
print("looking for an EEG stream...")
streams = resolve_stream('type', 'EEG')
inlet = StreamInlet(streams[0])
info = inlet.info()

# Retrieve channel names
ch = info.desc().child("channels").child("channel")
ch_names = []

num_channels = 2
sfreq = info.nominal_srate()

# Create a buffer for a 250 ms window
window_size = int(sfreq * 0.25)
buffer = Buffer(window_size, num_channels)

# Create a buffer for the first 30 seconds of median power values
median_power_buffer_size = bsl * 8  # 30 seconds, 8 windows per second
median_power_buffer = Buffer(median_power_buffer_size, num_channels)

# Create a second buffer to track the most recent 500 values outside the current processing window
history_buffer_size = 500
history_buffer = Buffer(history_buffer_size, num_channels)

num_new_samples_since_last_window = 0
window_count = 0
threshold_calculated = False
threshold = None
consecutive_below_threshold = np.zeros(num_channels)  # Array to track consecutive below-threshold windows

# LSL outlet for filtered EEG data
filtered_info = StreamInfo('FilteredEEG', 'EEG', num_channels, sfreq, 'float32', 'filtered_eeg123')
chns = filtered_info.desc().append_child("channels")
for chan_ix, label in enumerate(ch_names):
    ch = chns.append_child("channel")
    ch.append_child_value("label", label)
    ch.append_child_value("unit", "microvolts")
    ch.append_child_value("type", "EEG")
filtered_outlet = StreamOutlet(filtered_info)

# LSL outlet for power values
power_info = StreamInfo('Power', 'Power', num_channels, sfreq, 'float32', 'power123')
chns = power_info.desc().append_child("channels")
for chan_ix, label in enumerate(ch_names):
    ch = chns.append_child("channel")
    ch.append_child_value("label", label)
    ch.append_child_value("unit", "microvolts")
    ch.append_child_value("type", "EEG")
power_outlet = StreamOutlet(power_info)

# Create Marker stream info and outlet to send trigger when threshold is triggered
thresh_trig_info = StreamInfo('thresh_Markers', 'Markers', 1, 0, 'string', 'myuid_beh')
thresh_trig_outlet = StreamOutlet(thresh_trig_info)

# Find the index of channel 'C3'
ch_names = []
ch = info.desc().child("channels").child("channel")
for k in range(info.channel_count()):
    ch_names.append(ch.child_value("label"))
    ch = ch.next_sibling()
C3_index = ch_names.index('C3')

while True:
    # Get chunk
    chunk, timestamps = inlet.pull_chunk(timeout=0.0)
    
    if chunk:
        # Transpose the chunk to match the buffer shape (channels × samples)
        chunk = np.array(chunk).T
        buffer.add(chunk)
        history_buffer.add(chunk)
        num_new_samples_since_last_window += chunk.shape[1]
        
        # Check if it's time to process a new window
        if num_new_samples_since_last_window >= window_size // 2:
            # Get data for current window
            current_window = buffer.get_latest_window()
            
            # Extract left padding from the history buffer
            left_padding = history_buffer.get_latest_window()

            # Combine padding with the current window (only left padding needed for FIR)
            padded_signal = np.concatenate((left_padding, current_window), axis=1)
            
            # Bandpass filter
            filtered_data, avg_gd = bandpass_filter(padded_signal, 13, 20, sfreq, numtaps=101)
            
            # Remove padding after filtering
            filtered_data = filtered_data[:, history_buffer_size:]

            # Calculate which part of the window is new
            new_filtered_data = filtered_data[:, -num_new_samples_since_last_window:]
            # Stream only the new data
            filtered_outlet.push_chunk(new_filtered_data.T.tolist())
            
            # Calculate power for the processed window
            median_power = calculate_power(filtered_data, instantaneous=False)

            # Track median power values for the first 240 windows (30 seconds)
            if window_count < median_power_buffer_size:
                median_power_buffer.add(median_power[:, np.newaxis])
                window_count += 1

            elif not threshold_calculated:
                # Calculate the threshold as the 25th percentile of the median power values for each channel
                data = median_power_buffer.get_data()
                threshold = np.percentile(data, 25, axis=1)
                threshold_calculated = True
                print("Threshold calculated:", threshold)

            elif threshold_calculated:
                below_threshold = (median_power < threshold).astype(float)
                consecutive_below_threshold = (consecutive_below_threshold + below_threshold) * below_threshold
                trigger_channels = (consecutive_below_threshold >= 3).astype(float)

                # If channel C3 meets the trigger condition
                if consecutive_below_threshold[C3_index] >= 3:
                    thresh_trigger = "go"
                    thresh_trig_outlet.push_sample([thresh_trigger])

                    if vibFeedback:
                        toggle_vib('on')  # Turn the vibration motor on
                        time.sleep(.2)
                        toggle_vib('off')  # Turn the vibration motor off

                # Reset count for channels that did not meet the threshold
                consecutive_below_threshold *= below_threshold
            
            # Stream power data
            plotPower = np.repeat(median_power[:, np.newaxis], filtered_data.shape[1], axis=1)
            new_power_data = plotPower[:, -(window_size // 2):]  # Correctly identifying new data for power
            power_outlet.push_chunk(new_power_data.T.tolist())
            
            num_new_samples_since_last_window -= window_size // 2
