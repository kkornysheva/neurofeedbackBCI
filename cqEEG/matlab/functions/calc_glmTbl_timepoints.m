function [glmTblNew] = calc_glmTbl(what,settings,glmTbl,varargin)

switch(what)
    case 'beh'
        A = varargin{1};
        corrTrial= varargin{2};
        sub= varargin{3};
        freq= varargin{4};
        T = varargin{5};

        seqRTpointsA = [];
        seqIntpointsA = [];

        timingPrec = zeros(1, length(corrTrial));

        for i = 1:length(corrTrial)
            trials = corrTrial(i);
            if T.trialType(trials)==1 % Sequence trials
                timingPerformance=A.timing(trials,:)-A.cueTime(trials,:);
                timingPerformanceRT=A.timing(trials,1); %seq RT for 1st press
                timingPerformance=timingPerformance(2:4) ./ diff([0 A.cueTime(trials,2:4)]) * 100; % Deviation from target timing as percent of target interval
                timingPerformance=[0 timingPerformance];

                %%%% Points for sequence RT:
                if timingPerformanceRT>=0 & timingPerformanceRT<=200
                    seqRTpoints=5;
                elseif timingPerformanceRT>200 & timingPerformanceRT<=360
                    seqRTpoints=4;
                elseif timingPerformanceRT>360 & timingPerformanceRT<=480
                    seqRTpoints=3;
                elseif timingPerformanceRT>480 & timingPerformanceRT<=560
                    seqRTpoints=2;
                elseif timingPerformanceRT>560 & timingPerformanceRT<=600
                    seqRTpoints=1;
                elseif timingPerformanceRT>600
                    seqRTpoints=0;
                elseif timingPerformanceRT<0  % if they press a button prematurely before the digit cue (tZero)
                    seqRTpoints=0;
                elseif ~ismember(trials, corrTrial)
                    seqRTpoints=0;
                end;

                seqRTpointsA = [seqRTpointsA seqRTpoints];

                %%%% Points for deviation in percent of target interval:
                timingPerformanceInt=timingPerformance(2:4);
                if nanmean(abs(timingPerformanceInt))>=0 & nanmean(abs(timingPerformanceInt))<=10
                    seqIntpoints=5;
                elseif nanmean(abs(timingPerformanceInt))>10 & nanmean(abs(timingPerformanceInt))<=20
                    seqIntpoints=4;
                elseif nanmean(abs(timingPerformanceInt))>20 & nanmean(abs(timingPerformanceInt))<=30
                    seqIntpoints=3;
                elseif nanmean(abs(timingPerformanceInt))>30 & nanmean(abs(timingPerformanceInt))<=40
                    seqIntpoints=2;
                elseif nanmean(abs(timingPerformanceInt))>40 & nanmean(abs(timingPerformanceInt))<=50
                    seqIntpoints=1;
                elseif nanmean(abs(timingPerformanceInt))>50
                    seqIntpoints=0;
                elseif A.timing(trials,1)<0 % if they press a button prematurely before the Go cue
                    seqIntpoints=0;
                elseif ~ismember(trials, corrTrial)
                    seqIntpoints=0;
                end
                seqIntpointsA = [seqIntpointsA seqIntpoints];
                trialpoints(i)=seqIntpoints+seqRTpoints; % Log points
                timingPrec(i) = nanmean(abs(timingPerformanceInt));


            end
        end

        %seqIntpointsA contains the points for timing accuracy per trial
        %seqRTpointsA contains the points for RT per RT
        %trialpoints contains the overall points per trial and is equal to A.points

        % I'm outputting 'nanmean(abs(timingPerformanceInt))' to the glmTbl
        % as a more accurate measure of temporal precision per trial for now
        % - only for cond 1 & 2, not for thumb

        % Remove RT outliers
        RTs = A.probeRT(corrTrial);
        % outlier removal: 3 Standard Deviations
        meanRT = nanmean(RTs);
        stdRT = nanstd(RTs);
        outlierIdx = abs(RTs - meanRT) > (3 * stdRT);
        RTs(outlierIdx) = NaN;
%         numNotNaN = sum(~isnan(RTs));

        if settings.cond == 5
            tempTable = table(repelem(sub,length(corrTrial))',corrTrial,A.tempID(corrTrial),RTs,A.points(corrTrial),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),'VariableNames',{'subject','trial','condition','RT','points','timingPrec','power','pow0_300','pow400_700','pow800_1100','pow1200_1500','prodStart','velGo','velHiLo','maxERDfreq'});
        else
            tempTable = table(repelem(sub,length(corrTrial))',corrTrial,A.tempID(corrTrial),RTs,A.points(corrTrial),timingPrec',NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),NaN(length(corrTrial),1),'VariableNames',{'subject','trial','condition','RT','points','timingPrec','power','pow0_300','pow400_700','pow800_1100','pow1200_1500','prodStart','velGo','velHiLo','maxERDfreq'});
        end
% 
%         if settings.cond == 5
%             tempTable = table(repelem(sub, length(corrTrial))', corrTrial, A.tempID(corrTrial), RTs, A.points(corrTrial), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 'VariableNames', {'subject', 'trial', 'condition', 'RT', 'points', 'timingPrec', 't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 't10', 't11', 't12', 't13', 't14', 't15', 't16', 't17', 't18', 't19', 't20', 't21', 't22', 't23', 't24', 't25', 't26', 't27', 't28', 't29', 't30', 't31', 't32', 't33', 't34', 't35', 't36', 't37', 't38', 't39', 't40', 't41', 't42', 't43', 't44', 't45', 't46', 't47', 't48', 't49', 't50', 't51', 't52', 't53', 't54', 't55', 't56', 't57', 't58', 't59', 't60', 't61', 't62', 't63', 't64', 't65', 't66', 't67'});
%         else
%             tempTable = table(repelem(sub, length(corrTrial))', corrTrial, A.tempID(corrTrial), RTs, A.points(corrTrial), timingPrec', ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), NaN(length(corrTrial), 1), ...
%                 NaN(length(corrTrial), 1), ...
%                 'VariableNames', {'subject', 'trial', 'condition', 'RT', 'points', 'timingPrec', 't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 't10', 't11', 't12', 't13', 't14', 't15', 't16', 't17', 't18', 't19', 't20', 't21', 't22', 't23', 't24', 't25', 't26', 't27', 't28', 't29', 't30', 't31', 't32', 't33', 't34', 't35', 't36', 't37', 't38', 't39', 't40', 't41', 't42', 't43', 't44', 't45', 't46', 't47', 't48', 't49', 't50', 't51', 't52', 't53', 't54', 't55', 't56', 't57', 't58', 't59', 't60', 't61', 't62', 't63', 't64', 't65', 't66', 't67'});
%         end


        glmTblNew = [glmTbl.(freq);tempTable];

    case 'tfr'
        TFR = varargin{1};
        freq = varargin{2};
        velocity = varargin{3};

        if strcmp(freq,'maxERD')
            % Find frequency band with lowest avg power in preparation period
            meanPower = squeeze(nanmean(nanmean(nanmean(TFR.powspctrm(:, settings.chan, 3:10, 14:29),1),2),4));
            [minPower, minIx] = min(meanPower);
            maxERDfreqIx = 2 + minIx;
        end

        for i = 1:size(TFR.powspctrm,1)
            if strcmp(freq,'all')
                trialPower = nanmean(TFR.powspctrm(i,settings.chan,8:30,14:29),'all'); % trial i, selected channel(s), mu+beta range, preparation period
                trialPow0_300 = nanmean(TFR.powspctrm(i,settings.chan,8:30,14:17),'all'); % trial i, C3, all freqs, 0:300 ms
                trialPow400_700 = nanmean(TFR.powspctrm(i,settings.chan,8:30,18:21),'all');
                trialPow800_1100 = nanmean(TFR.powspctrm(i,settings.chan,8:30,22:25),'all');
                trialPow1200_1500 = nanmean(TFR.powspctrm(i,settings.chan,8:30,26:29),'all');
                trialProdStart = nanmean(TFR.powspctrm(i,settings.chan,8:30,30),'all');
            elseif strcmp(freq,'maxERD')
                trialPower = nanmean(TFR.powspctrm(i,settings.chan,maxERDfreqIx-2:maxERDfreqIx+2,14:29),'all'); % trial i, selected channel(s), mu+beta range, preparation period
                trialPow0_300 = nanmean(TFR.powspctrm(i,settings.chan,maxERDfreqIx-2:maxERDfreqIx+2,14:17),'all'); % trial i, C3, all freqs, 0:300 ms
                trialPow400_700 = nanmean(TFR.powspctrm(i,settings.chan,maxERDfreqIx-2:maxERDfreqIx+2,18:21),'all');
                trialPow800_1100 = nanmean(TFR.powspctrm(i,settings.chan,maxERDfreqIx-2:maxERDfreqIx+2,22:25),'all');
                trialPow1200_1500 = nanmean(TFR.powspctrm(i,settings.chan,maxERDfreqIx-2:maxERDfreqIx+2,26:29),'all');
                trialProdStart = nanmean(TFR.powspctrm(i,settings.chan,maxERDfreqIx-2:maxERDfreqIx+2,30),'all');
            else
                trialPower = nanmean(TFR.powspctrm(i,settings.chan,:,14:29),'all'); % trial i, C3, selected freq band, preparation period
                trialPow0_300 = nanmean(TFR.powspctrm(i,settings.chan,:,14:17),'all');
                trialPow400_700 = nanmean(TFR.powspctrm(i,settings.chan,:,18:21),'all');
                trialPow800_1100 = nanmean(TFR.powspctrm(i,settings.chan,:,22:25),'all');
                trialPow1200_1500 = nanmean(TFR.powspctrm(i,settings.chan,:,26:29),'all');
                trialProdStart = nanmean(TFR.powspctrm(i,settings.chan,:,30),'all');
            end
            glmTbl.(freq).power(find(isnan(glmTbl.(freq).power),1)) = trialPower;
            glmTbl.(freq).pow0_300(find(isnan(glmTbl.(freq).pow0_300),1)) = trialPow0_300;
            glmTbl.(freq).pow400_700(find(isnan(glmTbl.(freq).pow400_700),1)) = trialPow400_700;
            glmTbl.(freq).pow800_1100(find(isnan(glmTbl.(freq).pow800_1100),1)) = trialPow800_1100;
            glmTbl.(freq).pow1200_1500(find(isnan(glmTbl.(freq).pow1200_1500),1)) = trialPow1200_1500;
            glmTbl.(freq).prodStart(find(isnan(glmTbl.(freq).prodStart),1)) = trialProdStart;
            glmTbl.(freq).velGo(find(isnan(glmTbl.(freq).velGo),1)) = velocity.goCue(i);
            glmTbl.(freq).velHiLo(find(isnan(glmTbl.(freq).velHiLo),1)) = velocity.hiLo(i);

            if strcmp(freq,'maxERD')
                glmTbl.(freq).maxERDfreq(find(isnan(glmTbl.(freq).maxERDfreq),1)) = 10 + maxERDfreqIx;
            end
        end

        glmTblNew = glmTbl.(freq);

    case 'save'
        filepath = fullfile(settings.saveDirGLMvData);
        checkDir(filepath);
        save(fullfile(filepath, 'glmTbl.mat'),'glmTbl');

end
end