function [TFR] = bslCorr_TFR(settings,TFR)

    if ~settings.bslperTrial
        % 1) FieldTrip
        cfg = [];
        cfg.baseline     = settings.bsl_tRange;
        cfg.baselinetype = 'relchange';
        cfg.maskstyle    = 'saturation';
        TFR              = ft_freqbaseline(cfg, TFR);
        TFR.powspctrm    = TFR.powspctrm .* 100; % relchange -> percent change

    elseif settings.bslperTrial
%         % 2) Same as 1) but by hand (baseline correction for individual trials)
%         timeIx = TFR.time >= settings.bsl_tRange(1) & TFR.time <= settings.bsl_tRange(2);
%         for trial = 1:size(TFR.powspctrm, 1)
%             tfdata = TFR.powspctrm(trial,:,:,:);
%             siz = size(tfdata);
%             tfdata = reshape(tfdata,siz(2:end));
%             meanBSL = repmat(nanmean(tfdata(:,:,timeIx), 3), [1 1 size(tfdata, 3)]);
%             TFR.powspctrm(trial, :, :, :) =  ((tfdata - meanBSL) ./ meanBSL) .* 100;
%         end


%        % 3) Calculate average baseline over all trials of individual subject+condition and substract it from each trial
        timeIX = TFR.time >= settings.bsl_tRange(1) & TFR.time <= settings.bsl_tRange(2);
        mean_baseline = nanmean(nanmean(TFR.powspctrm(:, :, :, timeIX), 1), 4);
        mean_baseline_replicated = repmat(mean_baseline, [size(TFR.powspctrm, 1), 1, 1, size(TFR.powspctrm, 4)]);
        TFR.powspctrm = ((TFR.powspctrm - mean_baseline_replicated) ./ mean_baseline_replicated) .* 100;
    
%         % 4) Subtract the average baseline over all trials from the first 20 trials, 
%         % and from trial 21 onward subtract the average baseline over the previous 20 trials.
%         4.1) Subtract the average baseline over all trials from the first 20 trials
%         timeIX = TFR.time >= settings.bsl_tRange(1) & TFR.time <= settings.bsl_tRange(2);
%         mean_baseline_all_trials = nanmean(nanmean(TFR.powspctrm(:, :, :, timeIX), 4), 1);
%         mean_baseline_replicated_all_trials = repmat(mean_baseline_all_trials, [size(TFR.powspctrm, 1), 1, 1, size(TFR.powspctrm, 4)]); 
%         % Apply the percent change baseline correction to the first 20 trials
%         TFR.powspctrm(1:20, :, :, :) = ((TFR.powspctrm(1:20, :, :, :) - mean_baseline_replicated_all_trials(1:20, :, :, :)) ./ mean_baseline_replicated_all_trials(1:20, :, :, :)) .* 100;
%         % 4.2) Subtract the average baseline over the previous 20 trials from trial 21 onwards
%         for trial = 21:size(TFR.powspctrm, 1)
%             % Extract the baseline data for the current trial, channel, and frequency
%             baseline_data_previous_trials = nanmean(TFR.powspctrm((trial-20):(trial-1), :, :, timeIX), 4);
%             % Calculate the mean baseline across trials for the current trial
%             mean_baseline_previous_trials = nanmean(baseline_data_previous_trials, 1);
%             % Replicate the mean_baseline_previous_trials to match the size of the baseline_data
%             mean_baseline_replicated_previous_trials = repmat(mean_baseline_previous_trials, [1, 1, 1, size(TFR.powspctrm, 4)]);
%             % Apply the baseline correction to the current trial
%             TFR.powspctrm(trial, :, :, :) = ((TFR.powspctrm(trial, :, :, :) - mean_baseline_replicated_previous_trials) ./ mean_baseline_replicated_previous_trials) .* 100;
%         end
%     
%     
%         % Plot
%         histogram(TFR.powspctrm(:,:,:,30:45))
%         figure();plot(TFR.time,squeeze(nanmean(nanmean(TFR.powspctrm(:,111,:,:),1),3)))
%         figure();plot(TFR.time,squeeze(nanmean(nanmean(TFR.powspctrm(1,111,:,:),1),3)))
%         figure();plot(TFR.time,squeeze(nanmean(nanmean(TFR.powspctrm(20,111,:,:),1),3)))
%         figure();plot(TFR.time,squeeze(nanmean(nanmean(TFR.powspctrm(22,111,:,:),1),3)))
%         figure();plot(TFR.time,squeeze(nanmean(nanmean(TFR.powspctrm(51,111,:,:),1),3)))
%         figure();plot(TFR.time,squeeze(nanmean(nanmean(TFR.powspctrm(71,111,:,:),1),3)))
%     
%         % plot single trials and subject average
%         figure();
%         for i = 1:size(TFR.powspctrm,1)
%             plot(TFR.time,squeeze(TFR.powspctrm(i,111,6,:)))
%             hold on;
%         end
%         plot(TFR.time,nanmean(squeeze(TFR.powspctrm(:,111,6,:))),"Color","k","LineWidth",5)
%         meanBSL = nanmean(TFR.powspctrm(:,111,6,9:13),'all') % baseline mean



% Plot single-trial TFR ... why doesn't it work properly?
%         cfg              = [];
%         cfg.xlim         = [-0.5 4.5];
%         cfg.ylim         = [6 40];
% %                     cfg.zlim         = [-30 30];
%         cfg.channel      = 'all'; %C3-contralateral M1- Consider neighbours? - the default takes avg power across all channels
%         cfg.layout       = 'biosemi128.lay'; % Layout for the Biosemi high density EEG system @ Bangor
%         cfg.comment = 'no';   % Turn off any default commentary in the plot
% 
%         fig = figure();
%         ft_singleplotTFR(cfg, TFR);
% 
%         pause(0.001);
    end
end