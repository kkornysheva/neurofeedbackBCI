# Config
filepath = '..\\..\\rec\\simpleRT\\pilot\\sub-P001\\ses-S001\\eeg'
filename = 'sub-P001_ses-S001_task-simpleRT_run-001_eeg.xdf'
XDF_eegstreamname = "obci_eeg"
XDF_LSLmarkerstreamname = "lsl_markers"
photosensor_streamname = "photoresistor"

# Import
import os
import matplotlib.pyplot as plt
import numpy as np
import mne
from mnelab.io.xdf import read_raw_xdf
from pyxdf import match_streaminfos, resolve_streams
from scipy.signal import filtfilt, butter, hilbert
import time

# Load EEG data and markers from XDF file
def load_data_and_markers(filepath, filename, eegstreamname, photosensor_streamname, markerstreamname):
    fullfile = os.path.join(filepath, filename)
    streams = resolve_streams(fullfile)
    ix_eeg = match_streaminfos(streams, [{"name": eegstreamname}])[0]
    nominal_eeg_srate = [s["nominal_srate"] for s in streams if s["stream_id"]==ix_eeg][0]

    ix_photo = match_streaminfos(streams, [{"name": photosensor_streamname}])[0]

    # Create raw object
    raw = read_raw_xdf(fullfile, stream_ids=[ix_photo], prefix_markers=True,fs_new=nominal_eeg_srate)

    # Extract events
    events_lsl, event_id = mne.events_from_annotations(raw)

    return raw, events_lsl, event_id

# Load  EEG data and markers
raw, events, event_ids = load_data_and_markers(filepath, filename, XDF_eegstreamname, photosensor_streamname, XDF_LSLmarkerstreamname)

# Plot first 60 seconds of photosensor data
sfreq = raw.info['sfreq']
start_stop_seconds = np.array([0, 60])
start_sample, stop_sample = (start_stop_seconds * sfreq).astype(int)

# Extract photosensor channel by name
raw_selection = raw['photoresistor_0', start_sample:stop_sample]

# Plot raw data
t = raw_selection[1]
signal = raw_selection[0].T
plt.figure(1)
plt.plot(t, signal)
plt.title("raw photosensor data")
plt.xlabel("time [s]")
plt.ylabel("signal")

# Events
mne.viz.plot_events(events, raw.info['sfreq'])
# assert len(events[:,2]) == 2000

# plot photosensor and go cues
go_cue_event_id = event_ids['2-go']  # Adjust this based on the actual event ID for '1-go'
go_cue_events = events[events[:, 2] == go_cue_event_id] # Filter events to only keep the 'go' events
go_cue_times = go_cue_events[:, 0] / sfreq  # Sample index to time in seconds
raw_selection = raw['photoresistor_0'] # Plot photosensor data
t = raw_selection[1]
signal = raw_selection[0].T
plt.figure(1)
plt.plot(t, signal)
plt.title("raw photosensor data")
plt.xlabel("time [s]")
plt.ylabel("signal")
# Overlay vertical lines at 'go' cues
for go_time in go_cue_times:
    plt.axvline(x=go_time, color='red', linestyle='-', linewidth=0.5)  # thin red vertical line
plt.show()

# Check drift if I can figure out how to have hardware and software triggers simultaneoulsy
# ...

# Epoch
epochs = mne.Epochs(
    raw,
    events,
    event_id = go_cue_event_id,
    tmin = -0.005,
    tmax = .15,
    baseline=None,
    preload=True)

epochs.plot_image('photoresistor_0')

epochs_df = epochs.to_data_frame()

# Analyse timings
# Restrict data to 5-95 % of range between max and min for robustness against outliers
data_max = max(epochs_df['photoresistor_0']) * 0.95
data_min =  min(epochs_df['photoresistor_0']) + abs(max(epochs_df['photoresistor_0']) - min(epochs_df['photoresistor_0'])) * 0.05
range_minmax = abs(data_max-data_min)

# Thresholds at 10% and 75%
threshold10 = data_min + range_minmax * 0.1
threshold75 = data_min + range_minmax * 0.75

# Function to find supra-threshold indices
def supraThresh(epochs_df,threshold):
    supraThresh = [epochs_df['time'].iloc[ix] for ix, val in enumerate(epochs_df['photoresistor_0']) if val > threshold]
    return supraThresh

# RTs -------------------------------------------------
# Reaction time (= monitor input lag)
reactionTimes = [supraThresh(epochs_grouped,threshold10)[0] for ix, epochs_grouped in epochs_df.groupby(['epoch'])]
reactionTime = np.mean(reactionTimes)

# Reponse time
responseTimes = [supraThresh(epochs_grouped,threshold75)[0] for ix, epochs_grouped in epochs_df.groupby(['epoch'])]
responseTime = np.mean(responseTimes)

# Raise time
raiseTime = responseTime - reactionTime

# Stimulus duration --------------------------------------
# = first to last value above 75 % threshold
stimDurSamples = [supraThresh(epochs_grouped,threshold75) for ix, epochs_grouped in epochs_df.groupby(['epoch'])]
stimDur = [len(col) for col in stimDurSamples]
valsInEpochBeforeZero = epochs_df[epochs_df.time == 0].index[0]
stimStart = np.mean([col[0] for col in stimDurSamples]) # = response time
meanDuration = np.mean([col[-1] for col in stimDurSamples])
medianDuration = np.median([col[-1] for col in stimDurSamples])

# Plot
plt.figure(3)
epochs_df.set_index('time').groupby('epoch')['photoresistor_0'].plot(color="gray", linewidth=0.1)
epochs_df.set_index('epoch').groupby('time')['photoresistor_0'].agg(lambda x: np.mean(x.values.tolist(), axis=0)).plot(color="black", linewidth=2)
h_thresh10 = plt.axhline(threshold10, color="palegreen", linewidth=1.5, label='threshold 10%')
h_thresh75 = plt.axhline(threshold75, color="forestgreen", linewidth=1.5, label='threshold 75%')
plt.axvline(0, color="black", linewidth=1, label='trigger')
plt.axvline(reactionTime, color="black", linewidth=1)
plt.axvline(responseTime, color="black", linewidth=1)
# for i in np.arange(1,20):
#     plt.axvline(i*0.008, color="red", linewidth=.5)
h_reactiontime = plt.hlines(y=data_max, xmin=0, xmax=reactionTime, color="skyblue", label='reaction time: '+str(round(reactionTime*1000,2))+' ms')
h_raiseTime = plt.hlines(y=data_max-range_minmax*0.1, xmin=reactionTime, xmax=responseTime, color="aqua", label='raise time: '+str(round(raiseTime*1000,2))+' ms')
h_responsetime = plt.hlines(y=data_max-range_minmax*0.2, xmin=0, xmax=responseTime, color="blue", label='response time: '+str(round(responseTime*1000,2))+' ms')
h_stimDur = plt.hlines(y=data_max-range_minmax*0.235, xmin=stimStart, xmax=stimStart+meanDuration-responseTime, color="crimson", label='stimulus duration: '+str(round(meanDuration*1000,2))+' ms')
plt.title("switch black2white")
plt.xlabel("time [s]")
plt.ylabel("signal")
plt.legend(handles = [h_reactiontime, h_raiseTime, h_responsetime, h_stimDur, h_thresh75, h_thresh10], bbox_to_anchor=(.9,.9))
plt.show()

# Plot stimulus duration histogram
plt.figure(4)
plt.hist([col[-1] for col in stimDurSamples], bins=30)
plt.axvline(meanDuration, color="red", linewidth=0.8)
plt.axvline(medianDuration, color="yellow", linewidth=0.8)
plt.title("histogram of stimulus durations")
plt.xlabel("duration [s]")
plt.ylabel("number of stimuli")
plt.legend(["mean: "+str(round(meanDuration,2)),"median: "+str(round(medianDuration,2))], loc="upper left")
plt.show()

# Boxplot
# Compute distributions of reaction times, response times, and raise times across all epochs
reactionTimes = [supraThresh(epochs_grouped, threshold10)[0] for ix, epochs_grouped in epochs_df.groupby(['epoch'])]
responseTimes = [supraThresh(epochs_grouped, threshold75)[0] for ix, epochs_grouped in epochs_df.groupby(['epoch'])]
raiseTimes = [responseTimes[ix] - reactionTimes[ix] for ix in range(len(reactionTimes))]
box_data = [reactionTimes, responseTimes, raiseTimes]

plt.figure(figsize=(10, 6))
plt.boxplot(box_data, patch_artist=True, labels=['Reaction Time', 'Response Time', 'Raise Time'])
colors = ['lightblue', 'lightgreen', 'lightcoral']
for patch, color in zip(plt.gca().artists, colors):
    patch.set_facecolor(color)
plt.xlabel("Timing Metrics")
plt.ylabel("Time [s]")
plt.title("Distribution of Reaction, Response, and Raise Times")
plt.grid(False)  # Add a grid for better readability
plt.show()
