% permutation test

close all; clear; clc
addpath('./functions/');
addpath('../../toolboxes/permutest');

perm_prepOnly = 0;
perm_prodOnly = 0;

%% Setup
settings.saveDirTFRvData    = './data/eeg63/data'; % Folder to save tfr voltage data to
sbj                         = 1:18;
condID                      = [1,2,5]; % %slow==1, fast==2, thumb==5 from A.tempID - trigger value for thumb is 11
settings.freq               = {'theta','mu','lowBeta','highBeta','beta','maxERD','all'};
settings.subjectDir         = {'s05','s08','s11','s12','s13','s14','s16','s17','s18','s19','s20','s21','s22','s23','s24','s26','s27','s28'};
method                      = 'Morl';
settings.bslperTrial        = 1;
settings.bsl_tRange         = [-.5 0];

C3 = [111,113,113,108,110,110,115,113,114,112,115,108,115,109,111,108,108,113];
D16 = [108,110,110,105,107,107,112,110,111,109,112,105,112,106,108,105,105,110];
D26 = [118,120,120,115,117,117,122,127,121,119,122,115,122,128,118,114,115,120];
D28 = [120,122,122,117,119,119,124,120,123,121,124,117,124,115,120,116,117,122];

%% MAIN 
% load data
allTFR = cell(1,length(sbj));
for cond = condID(1:end)
    condition = getCondStr(cond);
    for freq = settings.freq
        freq = freq{:};
        Perm.(condition).(freq).(freq) = struct();
        for subIx = 1:length(sbj)
            sub = sbj(subIx);
            subTFR = fullfile(settings.saveDirTFRvData, 'sub/', settings.subjectDir{sub}, sprintf('TFR%s_%s_%s_%s.mat',method,freq,getCondStr(cond),settings.subjectDir{sub}));
            load(subTFR);
            % -------------------------------------------------------------
            % This can be run on eeg22 which isnt baseline corrected yet.
            % Results in the same thing as running eeg21 since we're bsl
            % correction on individual trials
%             TFR = rem_powOutliers(TFR);
%             TFR = bslCorr_TFR(settings,TFR);

            % This can be run on eeg31, which still is 4D to be run on Mike
            % X Cohen's code
            % adjust from singleTrial to subjectAverage
            TFR.powspctrm = squeeze(nanmean(TFR.powspctrm(:,:,:,:),1)); % mean across trials - chanxfreqxtime
            TFR.dimord = 'chan_freq_time'; % correct dimord according to the line above
            % -------------------------------------------------------------
            allTFR{sub} = TFR;
        end
        checkDir(fullfile(settings.saveDirTFRvData,'group/'))

        traceSub_all = [];
        bslSub_all = [];
        for i = 1:length(sbj)
            % calculate traceplot per subject
            if perm_prepOnly
                traceSub = squeeze(mean(mean(allTFR{1,i}.powspctrm(:,:,14:29),2),1)); % average over chan & freq dims
                times = allTFR{1,1}.time(14:29);
            elseif perm_prodOnly
                traceSub = squeeze(mean(mean(allTFR{1,i}.powspctrm(:,:,29:59),2),1)); % average over chan & freq dims
                times = allTFR{1,1}.time(29:59);
            else
                whichChan = [C3(i), D16(i), D26(i), D28(i)];
%                 traceSub = squeeze(mean(mean(allTFR{1,i}.powspctrm(:,:,14:59),2),1)); % average over chan & freq dims
%                 traceSub = squeeze(mean(allTFR{1,i}.powspctrm(whichChan(1),:,14:59),2)); % C3
                traceSub = squeeze(mean(mean(allTFR{1,i}.powspctrm(whichChan,:,14:59),2),1)); % Cluster4
                times = allTFR{1,1}.time(14:59);
            end

            bslSub = squeeze(mean(mean(allTFR{1,i}.powspctrm(:,:,9:14),2),1)); % average over chan & freq dims
            timesBsl = allTFR{1,1}.time(9:14);

            traceSub_all = [traceSub_all; traceSub'];
            bslSub_all = [bslSub_all; bslSub'];
        end

        % Equal options 1)
        traceGA = mean(traceSub_all,1);
        bslGA = mean(bslSub_all,1);
        % figure(); plot(TFR.time,traceGA);

        % Equal options 2)
%         cfg = [];
%         gaTFR = ft_freqgrandaverage(cfg,allTFR{:});
%         timeCourse = squeeze(mean(mean(gaTFR.powspctrm(:,:,:),2),1));
%         figure(); plot(gaTFR.time,timeCourse);
%         xlabel('time');
%         ylabel('power (percent change)');
%         title(sprintf('gaTFR traceplot %s %s',getCondStr(cond),freq))


        % Baseline --------------------------------------------------------
        % Option 1: Zeros
%         permBsl_fullSize = zeros(size(traceSub_all)); 
        % Option 2: Use real baseline (have to adapt size to be equal to data(=traceSub_all))
        numFullRepeats = floor(size(traceSub_all,2) / size(bslSub_all,2)); % Calculate full repeats
        remainder = mod(size(traceSub_all,2), size(bslSub_all,2)); % Calculate the remainder
        permBsl_fullSize = repmat(bslSub_all, 1, numFullRepeats); % Repeat baselineMatrix(=bslSub_all) for full repeats
        rng(0);
        if remainder > 0 % Randomly select columns from B to fill the remainder
            randomColumns = randi(size(bslSub_all,2), 1, remainder);
            permBsl_fullSize = [permBsl_fullSize, bslSub_all(:, randomColumns)];
        end



        % Cluster-based permutation test ----------------------------------

        % input traceplot per participant: subj x trace
        [Perm.(condition).(freq).clusters, Perm.(condition).(freq).p_values, Perm.(condition).(freq).t_sums, Perm.(condition).(freq).permutation_distribution] = ...
            permutest(traceSub_all',permBsl_fullSize','true',.05,1e4,'false');


        % Plot ------------------------------------------------------------

        figHandle = figure();
        hold on;


        % Option 1 --------------------------------------------------------
%         % Add shaded regions for significant clusters
%         for clusterIdx = 1:length(Perm.(condition).(freq).clusters)
%             pValue = Perm.(condition).(freq).p_values(clusterIdx);
%             
%             % Only display clusters with p-value <= 0.05
%             if pValue <= 0.05
%                 clusterTimes = times(Perm.(condition).(freq).clusters{clusterIdx});
%                 
%                 % Create vectors for fill function
%                 xValues = [clusterTimes, fliplr(clusterTimes)];
%                 yValues = [ones(size(clusterTimes))* 40, ones(size(clusterTimes))*30]; % Shading at the top
%                 
%                 % Draw shaded region for each cluster
%                 fill(xValues, yValues, 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%                 
%                 % Display rounded p-value as text
%                 roundedPValue = round(pValue, 3);
%                 text(mean(clusterTimes), 30, sprintf('p = %.3f', roundedPValue), ...
%                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
%             end
%         end
%         plot([0, 0], [-50, 50], 'k', 'LineWidth', 1); % Add vertical line at 0
%         plot([1.5, 1.5], [-50, 50], '--k', 'LineWidth', 1); % Add dotted vertical line at 1.5
%         hold off;
%         xlim([-.5, 4.5]);
%         ylim([-50, 50]);
%         xlabel('Time');
%         ylabel('Power (Percent Change)');
%         title(sprintf('GA TFR Traceplot %s %s', getCondStr(cond), freq));

        % Option 2 --------------------------------------------------------
        % Add bars and asterisks for significant clusters
        for clusterIdx = 1:length(Perm.(condition).(freq).clusters)
            pValue = Perm.(condition).(freq).p_values(clusterIdx);
            
            % Only display clusters with p-value <= 0.05
            if pValue <= 0.05
                clusterTimes = times(Perm.(condition).(freq).clusters{clusterIdx});
                
                % Create vectors for fill function
                xValues = [clusterTimes, fliplr(clusterTimes)];
                yValues = [ones(size(clusterTimes))*56, ones(size(clusterTimes))*55]; % Adjust the y position based on max value of traceGA
                
                % Draw shaded region for each cluster
                fill(xValues, yValues, 'k', 'FaceAlpha', 1, 'EdgeColor', 'none'); % 'k' for black bar
                
                % Determine the number of asterisks based on the p-value
                if pValue <= 0.001
                    significance = '***'; % p <= 0.001
                elseif pValue <= 0.01
                    significance = '**'; % p <= 0.01
                else
                    significance = '*'; % p <= 0.05
                end
                
                % Display asterisks at the top of each significant cluster
                text(mean(clusterTimes), 51, significance, ... % Adjust vertical position slightly above the bar
                     'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
            end
        end

        
        % Shaded error bars - std over sujects ----------------------------
        shadedErrorBar([timesBsl, times],[bslGA,traceGA],std([bslSub_all, traceSub_all]),'lineprops','b');

        
%         % Plot individual subjects
%         for i = 1:size(traceSub_all, 1)
%             plot([timesBsl,times], [bslSub_all(i, :),traceSub_all(i, :)], 'LineWidth', 1, 'Color', 'k');
%         end

        % Plot GA
        bslPlot = plot(timesBsl,bslGA,'LineWidth',4,'Color','b');
        curvePlot = plot(times,traceGA,'LineWidth',4,'Color','b');

        plot([0, 0], [-50, 60], 'k', 'LineWidth', 1); % Add vertical line at 0
        plot([1.5, 1.5], [-50, 60], '--k', 'LineWidth', 1); % Add dotted vertical line at 1.5
        plot([-.5, 4.5], [0, 0], ':k', 'LineWidth', 1); % Add dotted horizontal line at 0
        hold off;
        xlim([-.5, 4.5]);
        ylim([-50, 60]); % Adjust ylim based on maximum traceGA value
        xticks([0,1,2,3,4]); % Specify the x-ticks
        yticks([-40,-20,0,20,40,60]); % Specify the y-ticks
        xlabel('Time [s]');
        ylabel('Power [percent change]');
        title(sprintf('Permutation test %s %s', getCondStr(cond), freq));

%         exportgraphics(figHandle,sprintf('./figures/perm_%s_%s.eps',getCondStr(cond),freq),'ContentType','vector','BackgroundColor','none');
        saveas(gcf, fullfile('./figures', sprintf('perm_%s_%s.jpg',getCondStr(cond),freq)));

    end
end