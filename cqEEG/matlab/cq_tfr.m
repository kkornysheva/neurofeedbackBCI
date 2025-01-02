% - baseline correction
% - outlier removal
% - calculate TFRs
% - generate table for GLM (to load in other script)
% - plot figures
% - subject- and group-level analyses
% - ...

close all; clear; clc
dbstop if error
addpath('./functions/');

%% Initialize Fieldtrip
addpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/fieldtrip-20230422');
ft_defaults

%% settings parameters
settings.saveDat            = 1; % 1=yes, 0=no
settings.savePlt            = 1;
settings.saveGlmTbl         = 1;
settings.calcTFR            = 1; % 1 = yes, 0 = load saved TFR

settings.results            = {};

% settings.chan               = 111;
% settings.chan               = [109,111,119,120];
% settings.chan               = [96,104,106,108,111,113,115,118,120]; % D2, D10, D1sprintf("sub-P%03i/ses-S%03i/eeg",sub,sub)2, D14, D16, D19, D21, D26, D28

settings.bsl_tRange         = [-.5 0];
settings.bslperTrial        = 1; % bsl correction: 1 = on trials, 0 = on subjects

settings.theta.foi          = 4:7;
settings.mu.foi             = 8:12;
settings.lowBeta.foi        = 13:20;
settings.highBeta.foi       = 21:30;
settings.beta.foi           = 14:30;
settings.maxERD.foi         = 11:22; % look for maximal ERD per subject within alpha & lowBeta range, allow for boundaries
settings.all.foi            = 4:40;
settings.freq               = {'theta'};
settings.freq               = {'theta','mu','lowBeta','highBeta','beta','maxERD','all'};

settings.inDir              = '../../myrto/cqEEG/data/eeg/ICAcorrected'; % Folder containing final ICA-corrected EEG FT data files
settings.behDir             = '../../myrto/cqEEG/data/behaviour'; % Folder containing behavioural data output
settings.saveDirTFRvData    = './data/eeg63/data'; % Folder to save tfr voltage data to
settings.saveDirTFRvFigs    = './data/eeg63/figures'; % Folder to save tfr voltage figures to
settings.saveDirGLMvData    = './data/eeg63/data/glm'; % Folder to save glmTable to
settings.subjectDir         = {'s05','s08','s11','s12','s13','s14','s16','s17','s18','s19','s20','s21','s22','s23','s24','s26','s27','s28'};
settings.subjectnumber      = [5,8,11,12,13,14,16,17,18,19,20,21,22,23,24,26,27,28];

settings.fast.preparation   = [0 1.5];
settings.fast.production    = [1.5 3];
settings.slow.preparation   = [0 1.5];
settings.slow.production    = [1.5 4.5];
settings.thumb.preparation  = [0 1.5];
settings.thumb.production   = [1.5 2.5];
% NOTE: Myrto defined the production period as starting from the Go cue to 
% the end of a Sequence trial (t = 3 s) or of a Single press trial (t = 1 s).

settings.period             = {'preparation','production'};

settings.condID             = [1,2,5]; % %slow==1, fast==2, thumb==5 from A.tempID - trigger value for thumb is 11

settings.pointsFilter       = 1; % Points gained filter to identify 'correct' trials (normally 1, indicating spatially correct)

%% MAIN -------------------------------------------------------------------

% cq_tfr('subTFR',settings,1:18,'Morl','singleplot');
cq_tfr('subTFR',settings,1:18,'Morl','topoplot');
cq_tfr('gaTFR',settings,1:18,'Morl','singleplot');
cq_tfr('gaTFR',settings,1:18,'Morl','topoplot');
cq_tfr('gaTFR',settings,1:18,'Morl','traceplot');
cq_tfr('gaTFR',settings,1,'Morl','topoplotMarkers');

%% Function ---------------------------------------------------------------
function cq_tfr(what,settings,sbj,method,plottype)
% input options:
% what      'subTFR'    single subject TFR
%           'grTFR'     group TFR -> avg over all subjects
% settings       settingstings
% sbj       subject(s)
% method    'Hann'
%           'Morl'
% plottype  'singleplot'
%           'topoplot'
%           'traceplot'

switch(what)
    case 'subTFR'
        glmTbl = struct;
        for freq = settings.freq
            freq=freq{:};
            glmTbl.(freq) = table();
        end
        for subIx = 1:length(sbj)
            sub = sbj(subIx);
            eegFile   = fullfile(settings.inDir, settings.subjectDir{sub}, ['cdrfe_cqEEG_' settings.subjectDir{sub} '_ICA_corrected.mat']); %FT data
            load(eegFile);


            C3 = [111,113,113,108,110,110,115,113,114,112,115,108,115,109,111,108,108,113];
            D16 = [108,110,110,105,107,107,112,110,111,109,112,105,112,106,108,105,105,110];
            D26 = [118,120,120,115,117,117,122,127,121,119,122,115,122,128,118,114,115,120];
            D28 = [120,122,122,117,119,119,124,120,123,121,124,117,124,115,120,116,117,122];
            
            settings.chan = [C3(subIx), D16(subIx), D26(subIx), D28(subIx)];
%             settings.chan = [C3(subIx)];
            countRs = 1;

            for cond = settings.condID(1:end)
                settings.cond = cond;
                for freq = settings.freq
                    freq=freq{:};

                   % Load behavioral data
                    [corrTrial, glmTbl,pressT] = loadBeh(settings,sub,glmTbl,freq,length(data.trial));

                    if countRs
                        if sub == 3
                            pressT = pressT(1:end-1,1)
                        end
                        % move stimulus-locked to reponse-locked data
                        for i=1:length(pressT)
                            data.time{i} = data.time{i} - pressT(i,1)/data.fsample;
                        end
    %                     cfgResponse = [];
    %                     cfgResponse.offset = -pressT(:,1);
    %                     data = ft_redefinetrial(cfgResponse, data);
                        countRs = 0;
                    end

                    % TFR
                    if settings.calcTFR
                        % calculate TFR
                        TFR = calc_TFR_ft(what,settings,data,freq,corrTrial,method);

                        % Remove outliers
                        TFR = rem_powOutliers(TFR);

                        % Baseline correction on single trials
                        if settings.bslperTrial
                            TFR = bslCorr_TFR(settings,TFR);
                        end

                    else 
                        % load TFR from disk
                        TFR   = fullfile(settings.saveDirTFRvData, 'sub/', settings.subjectDir{sub}, sprintf('TFR%s_%s_%s_%s.mat',method,freq,getCondStr(settings.cond),settings.subjectDir{sub}));
                        load(TFR);
                    end

                    % Get channel with max ERD to look an individual differences
%                     [minVal,minIx] = min(nanmean(nanmean(nanmean(TFR.powspctrm(:,:,:,14:29),4),3),1));
%                     settings.results = [settings.results; {sub, settings.cond, freq, minVal, minIx, TFR.label(minIx)}];
                    

                    % calculate velocity of TFR
                    velocity = calc_TFRvel(TFR);

                    % glmTbl
                    glmTbl.(freq) = calc_glmTbl('tfr',settings,glmTbl,TFR,freq,velocity);
    
%                         % adjust from singleTrial to subjectAverage
%                         TFR.powspctrm = squeeze(nanmean(TFR.powspctrm(:,:,:,:),1)); % mean across trials - chanxfreqxtime
%                         TFR.dimord = 'chan_freq_time'; % correct dimord according to the line above
%                         
%                          % Baseline correction on subject average
%                         if ~settings.bslperTrial
%                             TFR = bslCorr_TFR(settings,TFR);
%                         end

                    % plot
                    plotTFR_ft(what,settings,TFR,freq,method,plottype,sub);
                    % save
                    saveTFR(what,settings,sub,method,TFR,freq);
                 end
            end

            % Raster plot
            plot_raster(glmTbl,pressT,freq, sub);
            close all;
        end
        if settings.saveGlmTbl 
            calc_glmTbl('save',settings,glmTbl);
        end

%         % Get channel with max ERD
%         maxERDchannels = cell2table(settings.results, 'VariableNames', {'Subject', 'Condition', 'Frequency', 'MinValue', 'MinIx', 'Label'});
%         filepath = fullfile(settings.saveDirGLMvData);
%         checkDir(filepath);
%         save(fullfile(filepath, 'erdChans.mat'),'maxERDchannels');
        
     case 'gaTFR'
        allTFR = cell(1,length(sbj));
        for cond = settings.condID(1:end)
            settings.cond = cond;
            for freq = settings.freq
                freq=freq{:};
                for subIx = 1:length(sbj)
                    sub = sbj(subIx);
                    subTFR = fullfile(settings.saveDirTFRvData, 'sub/', settings.subjectDir{sub}, sprintf('TFR%s_%s_%s_%s.mat',method,freq,getCondStr(settings.cond),settings.subjectDir{sub}));
                    load(subTFR);

                    % adjust from singleTrial to subjectAverage
                    TFR.powspctrm = squeeze(nanmean(TFR.powspctrm(:,:,:,:),1)); % mean across trials - chanxfreqxtime
                    TFR.dimord = 'chan_freq_time'; % correct dimord according to the line above

                    allTFR{sub} = TFR;
                end

                % calculate ga TFR
                cfg = [];
                if settings.saveDat
                    checkDir(fullfile(settings.saveDirTFRvData,'group/'))
                    cfg.outputfile = fullfile(settings.saveDirTFRvData,'group/', sprintf('gaTFR%s_%s_%s.mat',method,freq,getCondStr(settings.cond)));
                end
                gaTFR = ft_freqgrandaverage(cfg,allTFR{:});

                % This code can be run on eeg22 to achieve the same plots
                % as in cqEEG_manuscript_1.0 Figure 3 ---------------------
%                 cfg = [];
%                 cfg.baseline     = [-0.5 0];
%                 cfg.baselinetype = 'relchange';
%                 gaTFR = ft_freqbaseline(cfg, gaTFR);
%                 gaTFR.powspctrm    = gaTFR.powspctrm .* 100; % relchange -> percent change
                % ---------------------------------------------------------
% 
%                 if strcmp(plottype,'traceplot')
% 
% %                Options:
%                 % 1) plot ga timecourse averaged over the freq band
% %                 timeCourse = squeeze(nanmean(gaTFR.powspctrm(111,:,:),2));
% 
% 
% %               % 2) plot trace for each channel
%                 timeCourse = squeeze(nanmean(gaTFR.powspctrm(:,:,:),2));
% 
% %               % 3)plot mean over all channels (currently equal to
% %                 traceplot in plotTFR_ft()
% %                     timeCourse = squeeze(nanmean(nanmean(gaTFR.powspctrm(:,:,:),2),1));
% 
%                     figure();
%                     plot(gaTFR.time,timeCourse);
%                     xlabel('time');
%                     ylabel('power (percent change)');
%                     title(sprintf('gaTFR traceplot %s %s',getCondStr(settings.cond),freq))
%                 end
    
                plotTFR_ft(what,settings,gaTFR,freq,method,plottype);
            end
        end
    end
end


