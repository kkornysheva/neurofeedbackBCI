% Analysis of the simpleRT data:
% - light preprocessing
% - managing events
% - correcting for display delay (measured with photoresistor)
% - TFR
% - ...

close all; clear; clc;
cd('/rds/projects/k/kornyshk-kornyshevalab/martin/simpleRT/matlab');

% variables
sub = 1;

% initialize Fieldtrip ----------------------------------------------------
addpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/fieldtrip-20230422');
ft_defaults

% paths -------------------------------------------------------------------
bidsPath = sprintf("../data/raw/sub-P%03i/ses-S001/eeg",sub);
behEventsFile = fullfile(bidsPath,sprintf("../beh/sub-%03i_ses-001_task-simpleRT_run-001_events.tsv",sub));
behEventsFile = fullfile(bidsPath,sprintf("../beh/sub-%03i_ses-001_task-simpleRT_run-001_events.tsv",sub));
eegName = sprintf("sub-P%03i_ses-S001_task-simpleRT_run-001_eeg.xdf",sub);
if sub == 7
    eegName = sprintf("sub-P%03i_ses-S001_task-simpleRT_run-002_eeg.xdf",sub);
end
eventsName = sprintf("sub-%03i_ses-001_task-simpleRT_run-001_events.tsv",sub);
eegFile = fullfile(bidsPath,eegName);
eventsFile = fullfile(bidsPath,eventsName);
if ~isfile(eegFile)
    error('The specified EEG file does not exist: %s', eegFile);
elseif ~isfile(eventsFile)
    error('The specified events file does not exist: %s', eventsName);
end

% load data ---------------------------------------------------------------
[lslData, lslEvents] = xdf2fieldtrip(eegFile);
% split data for eeg ...
if isequal(lslData.label(4), {'1'}) % first three channels are analog
    eegRaw = lslData;  % Copy the original structure
    eegRaw.label = lslData.label(4:19);  % Select EEG labels
    eegRaw.trial = cellfun(@(x) x(4:19, :), lslData.trial, 'UniformOutput', false);  % Select EEG trials
    % ... and photoresistor
    photoRaw = lslData;
    photoRaw.label = lslData.label(1:3);
    photoRaw.trial = cellfun(@(x) x(1:3, :), lslData.trial, 'UniformOutput', false);
elseif isequal(lslData.label(4), {'4'}) % last three channels are analog
    eegRaw = lslData;  % Copy the original structure
    eegRaw.label = lslData.label(1:16);  % Select EEG labels
    eegRaw.trial = cellfun(@(x) x(1:16, :), lslData.trial, 'UniformOutput', false);
    photoRaw = lslData;
    photoRaw.label = lslData.label(17:19);
    photoRaw.trial = cellfun(@(x) x(17:19, :), lslData.trial, 'UniformOutput', false);
end

srate = lslData.hdr.Fs;  % Sampling rate in Hz

% Add channel labels
eegRaw.label = {'O1', 'CP3', 'CP5', 'C1', 'C3', 'C5', 'CP1', 'Fp1', ...
              'O2', 'CP4', 'CP6', 'C2', 'C4', 'C6', 'CP2', 'Fp2'}';

%% preproc eeg 

% plot some raw eeg data 
figure; plot(eegRaw.time{1,1}(1:100000),eegRaw.trial{1,1}(1,1:100000)); 

cfgEEG = [];
cfgEEG.reref = 'yes';
cfgEEG.refchannel = 'all';  % common average reference
cfgEEG.hpfilter = 'yes';
cfgEEG.hpfreq = 0.5;  % Cut-off frequency for highpass filter
eegPreproc = ft_preprocessing(cfgEEG, eegRaw); % Apply the preprocessing to the data loaded with xdf2fieldtrip

% confirm hp filter with a little plot 
figure; plot(eegPreproc.time{1,1}(1:100000),eegPreproc.trial{1,1}(1,1:100000)); 

%% visualise photosensor

% discard unused analog channels 
cfgPhoto = [];
cfgPhoto.channel = 1; 
photoRaw = ft_preprocessing(cfgPhoto, photoRaw);

figure;
plot(photoRaw.time{1,1},photoRaw.trial{1,1}(1,:));
hold on;
goCueIndices = find(strcmp({lslEvents.value}, 'go'));  % Find indices of 'go' events
for i = 1:length(goCueIndices)
    goTime = lslEvents(goCueIndices(i)).timestamp; % Get the timestamp of each 'go' cue
    line([goTime goTime], ylim, 'Color', 'r', 'LineWidth', 0.5);  % Overlay red line at go cue time
end
hold off;

%% Events

% Update the FieldTrip events with sample indices calculated from timestamps
% (relative to the first sample !!!)
for i = 1:length(lslEvents)
    lslEvents(i).sample = round((lslEvents(i).timestamp - eegPreproc.hdr.FirstTimeStamp) * srate);  % Convert timestamp to sample index
end

% Read events from file
evt = ft_read_tsv(eventsFile);
% Extract events from block 1, which are for the first block (training, only beh, no eeg)
events_beh = evt(evt.block == 0 | evt.block == 1, :);
% Save the table 'events_beh' as a .tsv file
writetable(events_beh, behEventsFile, 'FileType', 'text', 'Delimiter', '\t');
% Extract events from block 1, which are for the first block (training, only beh, no eeg)
events_eeg = evt(evt.block == 0 | evt.block > 1, :);
% Save the table 'events_beh' as a .tsv file
writetable(events_eeg, eventsFile, 'FileType', 'text', 'Delimiter', '\t');
% Make compatible with FieldTrip expectations
events_eeg.Properties.VariableNames{'event'} = 'value';
events_eeg = table2struct(events_eeg);
[events_eeg.type] = deal('marker');
[events_eeg.offset] = deal([]);
sample_values = num2cell(round([events_eeg.onset] * srate));
[events_eeg.sample] = deal(sample_values{:});
events_eeg(1) = [];

if sub == 2 || sub == 4
    lslEvents(end) = [];
elseif sub == 7
    % Remove events from block 1
    events_eeg = events_eeg(232:1800);
    lslEvents = lslEvents(1:1569);
end

% confirm that the difference between events and events_eeg is minimal
y1 = [events_eeg.onset]; % in seconds
y2 = [lslEvents.timestamp]; % in seconds
figure; plot(diff(y1)-diff(y2))

% since these are identical, we can copy over the timestamps from events to events_eeg
[events_eeg.sample] = deal(lslEvents.sample);

% Create final events dataframe
events = events_eeg;

% Indicate whether a trial is correct or not for all rows of each trial
for i = 1:length(events)
    if ~isnan(events(i).corr_trial) % check if corr_trial is 0 or 1 (i.e., skip if NaN)
        trial_value = events(i).trial;
        block_value = events(i).block;
        
        % Find all entries with the same trial and block combination
        match_idx = find([events.trial] == trial_value & [events.block] == block_value);
        
        % Copy the corr_trial value to all matching entries
        for idx = match_idx
            events(idx).corr_trial = events(i).corr_trial;
        end
    end
end

% Re-calculate response times based on onset column
% (this is more accurate than the ones that are in the reponse_time column!!!)
[events.rt] = deal(NaN); % initialize rt field with NaNs

for b = unique([events.block])
    for t = unique([events.trial])
        % Find the 'go' event for this block and trial
        go_idx = find(strcmp({events.value}, 'go') & [events.block] == b & [events.trial] == t);
    
        if ~isempty(go_idx)
            go_onset = events(go_idx).onset;  % Onset of the 'go' event
            
            % Find all 'r', 'w', 'e', or 'a' events for this block and trial
            rwea_idx = find(ismember({events.value}, {'r', 'w', 'e', 'a'}) & [events.block] == b & [events.trial] == t);
            
            % For each rwea event, calculate the reaction time
            for i = rwea_idx
                events(i).rt = events(i).onset - go_onset;  % Subtract go onset from rwea onset
            end
        end
    end
end

% Remove incorrect trials
events = events([events.corr_trial] == 1);

% Remove all rows with RT's > 3
valid_rows = [events.response_time] <= 3  | isnan([events.response_time]);
events = events(valid_rows);

length(find(strcmp({events.value}, 'go')))

%% epoching

cfgEvts = [];
cfgEvts.dataset = '';
cfgEvts.headerfile  = '';
cfgEvts.hdr = eegPreproc.hdr;
cfgEvts.event = events;
cfgEvts.trialdef.eventtype = 'marker';
cfgEvts.trialdef.eventvalue = 'go'; % 'go' for stimulus-locked and 'r' for response-locked analysis
cfgEvts.trialdef.prestim = 1;
cfgEvts.trialdef.poststim = 3.5;
cfgEvts.trialfun = 'ft_trialfun_general';
evts = ft_definetrial(cfgEvts);
data_epoch = ft_redefinetrial(evts, eegPreproc);

cfgEvts = [];
cfgEvts.dataset = '';
cfgEvts.headerfile  = '';
cfgEvts.hdr = eegPreproc.hdr;
cfgEvts.event = events;
cfgEvts.trialdef.eventtype = 'marker';
cfgEvts.trialdef.eventvalue = 'go';
cfgEvts.trialdef.prestim = 1;
cfgEvts.trialdef.poststim = 3.5;
cfgEvts.trialfun = 'ft_trialfun_general';
evts = ft_definetrial(cfgEvts);
photo_epoch = ft_redefinetrial(evts, photoRaw);

% also epoch around the 30 s baseline interval so I can use this for bsl
% correction
cfgEvts = [];
cfgEvts.dataset = '';
cfgEvts.headerfile  = '';
cfgEvts.hdr = eegPreproc.hdr;
cfgEvts.event = events_eeg;
cfgEvts.trialdef.eventtype = 'marker';
cfgEvts.trialdef.eventvalue = 'bslStart';
cfgEvts.trialdef.prestim = -5;
cfgEvts.trialdef.poststim = 35;
cfgEvts.trialfun = 'ft_trialfun_general';
evtsBsl = ft_definetrial(cfgEvts);

data_bsl = ft_redefinetrial(evtsBsl, eegPreproc);

%% Remove monitor RTs from eeg and photo data

% Plot single trial photo data
figure;
for i = 1:length(photo_epoch.trial)
    plot(photo_epoch.time{1}(126:139),photo_epoch.trial{i}(126:139));
    hold on;
end
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('unaligned single trial photosensor data');

% -------------------------------------------------------------------------
% Easy solution but retains 8 ms uncertainty around samples caused by 
% 125 Hz srate.

% % Correct for display latency (and other minor things that may contribute)
% % Get vector of monitor RTs
% monitorRTs = [];
% for i = 1:length(photo_epoch.trial)
%     trial_data = photo_epoch.trial{i};
%     idx = find(trial_data > 700, 1);  % Find index of the first value in the trial that is greater than 700
%     monitorRTs(i) = photo_epoch.time{i}(idx);  % Use the time vector of the current trial
% end
% 
% % Move data to adjust for monitor RT
% for i = 1:length(monitorRTs)
%     photo_epoch.time{i} = photo_epoch.time{i} - monitorRTs(i) + 0.004;  % Subtract the RT from each time vector
%     % we have 8ms uncertainty due to the 125 Hz srate. So I center the data on 0 +- 4 ms.
% 
%     % Do the same for our eeg data!
%     data_epoch.time{i} = data_epoch.time{i} - monitorRTs(i) + 0.004;
% end
% 
% % Plot the aligned single trial photo data
% figure;
% for i = 1:length(photo_epoch.trial)
%     plot(photo_epoch.time{i}(126:139), photo_epoch.trial{i}(126:139));
%     hold on;
% end
% hold off;
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('single trial photosensor data (medium precision)');
% % -> Now the data at an amplitude of 700 should all be between -4 to 4 ms!

% -------------------------------------------------------------------------
% How could I make this even more accurate? I.e. how could I remove the 8
% ms uncertainty around samples? To align all trials absolutely perfectly
% at my threshold value (700 - this is arbitrary - could e.g. go for 
% 75th percentile between central 90% of probability mass) I can't look for 
% the timing of the first sample that is larger than the threshold but I 
% would probably have to fit a function on each trial to have a continuous 
% function (where I can determine the x for every arbitrary y value) from 
% which I can continue the exact timepoint where the threshold is crossed.

% Improved accuracy by interpolating around the threshold crossing

% Initialize monitorRTs
monitorRTs = [];

% Loop through each trial in photo_epoch
for i = 1:length(photo_epoch.trial)
    trial_data = photo_epoch.trial{i};
    time_data = photo_epoch.time{i};
    
    % Find the first index where the value exceeds 700
    idx = find(trial_data > 700, 1);  % Find index of the first value greater than 700
    
    if ~isempty(idx) && idx > 1
        % Get two points: one before and one after the threshold crossing
        x1 = time_data(idx-1);
        x2 = time_data(idx);
        y1 = trial_data(idx-1);
        y2 = trial_data(idx);

        % Solve for the exact time when y = 700 (threshold crossing)
        % OPTION 1: Linear interpolation (analytic solution)
        slope = (y2 - y1) / (x2 - x1);
        intercept = y1 - slope * x1;
        exact_time = (700 - intercept) / slope;
        % OPTION 2: Spline Fitting (numerical solution)
        exact_time = interp1([y1, y2], [x1, x2], 700, 'spline');
        % -> both options work well here - practically doesn't make a
        % difference
        
        % Store the exact monitor RT for this trial
        monitorRTs(i) = exact_time;
    else
        monitorRTs(i) = NaN;  % In case no crossing is found
    end
end

% Move data to adjust for monitor RT with high accuracy
for i = 1:length(monitorRTs)
    if ~isnan(monitorRTs(i))
        % Adjust both photo_epoch and data_epoch time vectors
        photo_epoch.time{i} = photo_epoch.time{i} - monitorRTs(i);
        data_epoch.time{i} = data_epoch.time{i} - monitorRTs(i);
    end
end

% Plot the aligned single trial photo data
figure;
for i = 1:length(photo_epoch.trial)
    plot(photo_epoch.time{i}(126:139), photo_epoch.trial{i}(126:139));
    hold on;
end
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('aligned single trial photosensor data (high precision)');


%% Subtract monitor RTs from behavioral RTs

monitor_idx = 1;  % Index to track monitorRTs

for i = 1:length(events)
    if ismember(events(i).value, {'r', 'w', 'e', 'a'})
        % Subtract the corresponding monitor RT from the event's rt value
        events(i).rt = events(i).rt - monitorRTs(monitor_idx);  
        
        % Check if this is the last event in the set ('a' indicates end of one group)
        if strcmp(events(i).value, 'a')
            monitor_idx = monitor_idx + 1;  % Move to the next monitorRT
        end
    end
end

%% Inspect behavioral data

% Number of correct trials
disp(sprintf('number of correct trials: %i',sum(strcmp({events.value},'go'))));

% Response time
r_indices = strcmp({events.value},'r');
r_rt_values = [events(r_indices).rt];
% mean, median, std
mean_rt = mean(r_rt_values);
median_rt = median(r_rt_values);
std_r_rt = std(r_rt_values);
% histogram
figure;
histogram(r_rt_values, 30);
hold on;
xline(mean_rt, 'r:', 'LineWidth', 2); 
xline(median_rt, 'g:', 'LineWidth', 2);
xline(mean_rt + std_r_rt, 'k:', 'LineWidth', 1.5);  % Mean + 1*std (blue dotted line)
xline(mean_rt - std_r_rt, 'k:', 'LineWidth', 1.5);  % Mean - 1*std (blue dotted line)
xlabel('reaction time (s)');
ylabel('n');
title('reaction time');
legend({'reaction times', sprintf('mean: %.3f',mean_rt), sprintf('median: %.3f',median_rt), sprintf('std: %.3f', std_r_rt)});
hold off;

% Movement time
a_indices = strcmp({events.value},'a');
a_rt_values = [events(a_indices).rt];
movement_times = [];
for i = 1:length(a_rt_values)
    movement_times(i) = a_rt_values(i) - r_rt_values(i);
end
% mean, median, std
mean_mt = mean(movement_times);
median_mt = median(movement_times);
std_mt = std(movement_times);
% histogram
figure;
histogram(movement_times, 30);
hold on;
xline(mean_mt, 'r:', 'LineWidth', 2); 
xline(median_mt, 'g:', 'LineWidth', 2);
xline(mean_mt + std_mt, 'k:', 'LineWidth', 1.5);  % Mean + 1*std (blue dotted line)
xline(mean_mt - std_mt, 'k:', 'LineWidth', 1.5);  % Mean - 1*std (blue dotted line)
xlabel('movement time (s)');
ylabel('n');
title('movement time');
legend({'movement times', sprintf('mean: %.3f',mean_mt), sprintf('median: %.3f',median_mt), sprintf('std: %.3f', std_mt)});
hold off;

%% Plot beta ERP

% cfgEEG = [];
% cfgEEG.reref = 'yes';
% cfgEEG.refchannel = 'all';  % Reference to the average of all channels
% cfgEEG.bpfilter = 'yes';
% cfgEEG.bpfreq = [13 20];  % Cut-off frequencies for bandpass filter
% eegPreprocBP = ft_preprocessing(cfgEEG, eegRaw); % Apply the preprocessing to the data loaded with xdf2fieldtrip
% data_epoch = ft_redefinetrial(evts, eegPreprocBP);
% 
% % Perform timelock analysis to compute the ERP
% cfgERP = [];
% cfgERP.keeptrials = 'no';  % Average across trials
% % cfgERP.channel = 13;
% timelock = ft_timelockanalysis(cfgERP, data_epoch);
% 
% % Plot the ERP
% cfgPlotERP = [];
% cfgPlotERP.xlim = [-0.5 3];  % Time range for plotting (matching prestim and poststim)
% cfgPlotERP.ylim = 'maxmin';    % Auto-scale y-axis
% cfgPlotERP.channel = 13;
% 
% ft_singleplotER(cfgPlotERP, timelock);

%% Plot single trial EEG data

figure;
for i = 1:length(data_epoch.trial)
    plot(data_epoch.time{1},data_epoch.trial{i});
    hold on;
end
hold off;

%% Plot periodogram

% Configure the frequency analysis for Welch's method
cfg = [];
cfg.method     = 'mtmfft';    % Use FFT-based method for frequency analysis
cfg.output     = 'pow';       % Compute power spectrum
cfg.taper      = 'hanning';   % Use Hanning taper for Welch's method
cfg.foilim     = [1 40];      % Frequency range of interest
cfg.pad        = 'maxperlen'; % Pad to the maximum segment length
cfg.keeptrials = 'no';        % Average across trials
cfg.channel    = 'all';       % Use all available channels (can specify 'EEG' or other channels if needed)

freq = ft_freqanalysis(cfg, data_epoch);

% Plot the Welch periodogram using ft_singleplotER
cfg = [];
cfg.channel = 'all';  % Plot all channels (change to specific channels as needed)
cfg.xlim    = [1 40]; % X-axis (frequency) range
cfg.ylim    = 'maxabs'; % Adjust y-axis to the maximum absolute value
figure;
ft_singleplotER(cfg, freq);
title('Welch Periodogram');
xlabel('Frequency (Hz)');
ylabel('Power');


%% Reject bad trials

cfg        = [];
cfg.metric = 'zvalue';
cfg.method = 'summary';
% data_rejbad = ft_rejectvisual(cfg,data_epoch);

% do it automatically after doing it once via the gui above

if sub == 1
    bad_trials = [77,86,100,158,168,176];
elseif sub == 2
    bad_trials = [11,68,69,153];
elseif sub == 5
    bad_trials = [17,57,59,216];
else
    bad_trials = [];
end
cfg = [];
cfg.trials = setdiff(1:length(data_epoch.trial), bad_trials);
data_rejbad = ft_redefinetrial(cfg, data_epoch);

% Plot single trial data AGAIN - have a look at how the trial removal
% changes our data
% figure;
% for i = 1:length(data_rejbad.trial)
%     plot(data_rejbad.time{1},data_rejbad.trial{i});
%     hold on;
% end
% hold off;


% doing it by hand, but doesn't work as well ------------------------------
% % Initialize parameters
% time_threshold = 1.5;  % Time after which to search for outliers
% z_threshold = 3;       % Z-score threshold for detecting outliers
% 
% % Find the index corresponding to the time threshold
% time_index = find(data_epoch.time{1} >= time_threshold, 1);
% 
% % Concatenate data across trials for the time points after the threshold
% data_after_threshold = cellfun(@(x) x(:, time_index:end), data_epoch.trial, 'UniformOutput', false);
% data_matrix = cat(3, data_after_threshold{:});  % Convert to 3D matrix (channels x time x trials)
% 
% % Flatten the data across channels and time for each trial
% flattened_trials = squeeze(max(abs(data_matrix), [], [1, 2]));  % Find the maximum absolute value in each trial after 1.5s
% 
% % Compute z-scores of the maximum values across trials
% z_scores = (flattened_trials - mean(flattened_trials)) / std(flattened_trials);
% 
% % Find trials with z-scores beyond the threshold
% outlier_trials = find(abs(z_scores) > z_threshold);
% 
% % Display outlier trials
% fprintf('Outlier trials detected: %s\n', strjoin(arrayfun(@num2str, outlier_trials, 'UniformOutput', false), ', '));
% 
% % Plot the data with outliers highlighted
% figure;
% for i = 1:length(data_epoch.trial)
%     if ismember(i, outlier_trials)
%         plot(data_epoch.time{1}, data_epoch.trial{i}, 'r');  % Plot outliers in red
%     else
%         plot(data_epoch.time{1}, data_epoch.trial{i}, 'b');  % Plot normal trials in blue
%     end
%     hold on;
% end
% hold off;
% title('Detected Outlier Trials (in Red)');
% xlabel('Time (s)');
% ylabel('Amplitude (\muV)');
% 
% % Remove the outlier trials from data_epoch
% good_trials = setdiff(1:length(data_epoch.trial), outlier_trials);  % Find non-outlier trials
% 
% % Create a new data_epoch structure with only good trials
% data_rejbad = ft_redefinetrial(struct('trials', good_trials), data_epoch);

%% ICA
% cfg = [];
% cfg.method = 'runica'; % EEGLAB implementation
% comp = ft_componentanalysis(cfg, data_rejbad); % ICA decomposition
% 
% % plot the components for visual inspection
% figure
% cfg = [];
% cfg.component = 1:16;       % specify the component(s) that should be plotted
% cfg.layout    = 'obci.lay'; % specify the layout file that should be used for plotting
% cfg.comment   = 'no';
% ft_topoplotIC(cfg, comp)
% 
% % plot the timecourse
% cfg = [];
% cfg.layout = 'obci.lay'; % specify the layout file that should be used for plotting
% cfg.viewmode = 'component';
% ft_databrowser(cfg, comp)
% 
% % remove the bad components and backproject the data
% cfg = [];
% cfg.component = [11,12,15]; % to be removed component(s)
% data_rejbad = ft_rejectcomponent(cfg, comp, data_rejbad)
% ........................................................................



%% Automatic artifact rejection
% % EOG
% cfg = [];
% 
% % channel selection, cutoff and padding
% cfg.artfctdef.zvalue.channel     = data_rejbad.label;
% cfg.artfctdef.zvalue.cutoff      = 4;
% cfg.artfctdef.zvalue.trlpadding  = 0;
% cfg.artfctdef.zvalue.artpadding  = 0.1;
% cfg.artfctdef.zvalue.fltpadding  = 0;
% 
% % algorithmic parameters
% cfg.artfctdef.zvalue.bpfilter   = 'yes';
% cfg.artfctdef.zvalue.bpfilttype = 'but';
% cfg.artfctdef.zvalue.bpfreq     = [2 15];
% cfg.artfctdef.zvalue.bpfiltord  = 4;
% cfg.artfctdef.zvalue.hilbert    = 'yes';
% 
% % feedback
% cfg.artfctdef.zvalue.interactive = 'yes';
% 
% [cfg, artifact_eog] = ft_artifact_zvalue(cfg,data_rejbad);
% 
% 





%% TFR: Morlet wavelet 

cfgTFR              = [];
cfgTFR.output       = 'pow';
cfgTFR.channel      = 'all';
cfgTFR.trials       = 'all'; % ft_freqanalysis will read in all trials from level data but will analyse only trlNr
cfgTFR.method       = 'wavelet';
cfgTFR.width        = 7; 
cfgTFR.toi          = -0.5:.1:3; % in sliding windows of 0.1 s     
cfgTFR.foi          = 1:40;
cfgTFR.keeptrials   = 'yes';
cfgTFR.t_ftimwin    = ones(length(cfgTFR.foi),1).*0.1; % length of time window = 0.1 sec
TFR                 = ft_freqanalysis(cfgTFR, data_rejbad); % TFR.powspctrm contains the power for chanxfreqxtime

%% Bsl corr

% OPTION 1: standard per trial
cfgBsl = [];
cfgBsl.baseline     = [-.5 0];
cfgBsl.baselinetype = 'relchange';
cfgBsl.maskstyle    = 'saturation';
TFR                 = ft_freqbaseline(cfgBsl, TFR);
TFR.powspctrm       = TFR.powspctrm .* 100; % relchange -> percent change

% OPTION 2: use 30s rest period(s) as baseline
% Bsl corr with 30 s rest period instead of single trial bsl
% extract data from a rest period
% TFR on bsl data
% cfg = [];
% cfg.output     = 'pow';
% cfg.channel    = 'all';
% cfg.trials     = [1,2,3];
% cfg.method     = 'wavelet';
% cfg.width      = 7; 
% cfg.toi        = 1:0.1:30;
% cfg.foi        = 1:40;
% cfg.keeptrials = 'no';
% cfg.t_ftimwin = ones(length(cfg.foi),1).*0.1; % length of time window = 0.1 sec
% TFR_bsl = ft_freqanalysis(cfg, data_bsl);
% % Compute the mean baseline power across the entire baseline period
% mean_bsl_power = nanmean(TFR_bsl.powspctrm, 3);
% % Manually apply the relative change correction
% for chan = 1:size(TFR.powspctrm, 2)  % Loop over channels
%     for freq = 1:size(TFR.powspctrm, 3)  % Loop over frequencies
%         baseline_value = mean_bsl_power(chan, freq);  % Mean power for this channel and frequency
%         TFR.powspctrm(:, chan, freq, :) = ...
%             (TFR.powspctrm(:, chan, freq, :) - baseline_value) / baseline_value;
%     end
% end
% TFR.powspctrm = TFR.powspctrm * 100;  % Convert to percentage change

%% plot TFR
cfgPlotTFR = [];
cfgPlotTFR.xlim         = [-0.5 3];
cfgPlotTFR.ylim         = [2 30];
cfgPlotTFR.zlim         = [-50 150];
cfgPlotTFR.channel      = 13; % C4

fig = figure();
ft_singleplotTFR(cfgPlotTFR, TFR);
hold on;
plot([0, 0], [2, 30], 'k', 'LineWidth', 1); % Add vertical line at 0
% plot([1.5, 1.5], [6, 40], '--k', 'LineWidth', 1); % Add dotted vertical line at 1.5
hold off;

%% plot "beta power erp" (and individual trials)

% average across channels and frequencies:
individual_trials = squeeze(nanmean(TFR.powspctrm(:, 13, 13:20, :),3));  % Trials x Time over low beta

% plot individual trials
figure;
hold on;
for i = 1:size(individual_trials, 1)
    plot(TFR.time, individual_trials(i, :), 'Color', [0.8, 0.8, 0.8]);  % light gray
end

% average across trials
betaline = mean(individual_trials, 1); 
plot(TFR.time, betaline, 'k-', 'LineWidth', 2);
xlim([-.5 3])
% ylim([-40 150])
xlabel('Time (s)');
ylabel('Power');
title('low beta power');
hold off;


%% topoplot

% !!! I added the custom layout of our headset myself to be able to use it
% here! How? I copied the 10-5 layout in
% 'toolboxes\fieldtrip-20230422\template\layout', removed all the channels
% that we aren't using, and saved it as obci.lay in the same folder. !!!

cfg              = [];
cfg.xlim         = [0.3 1];
% cfg.zlim         = [-Inf 0.1];
cfg.colorbar     = 'yes';
cfg.marker       = 'on';
cfg.layout       = 'obci.lay';
cfg.comment      = 'no';   % Turn off any default commentary in the plot

fig = figure();
ft_topoplotTFR(cfg, TFR);


%% ...



