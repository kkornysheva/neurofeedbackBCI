function [corrTrial,glmTbl,pressT] = loadBeh(settings,sub,glmTbl,freq,dataLen)

    dataBeh          = fullfile(settings.behDir, settings.subjectDir{sub}, 'behDataEEG.mat');
    dataTarget       = fullfile(settings.behDir, settings.subjectDir{sub}, 'targetEEG.mat');
    load(dataBeh);          % Load 'A' structure containing the behavioural data
    load(dataTarget);       % Load 'T' structure containing more behavioural data

    % Remove respective bad trials from behavioural data structs:
    badTrials        = fullfile(settings.inDir, settings.subjectDir{sub}, 'badTrials.mat');
    load(badTrials);
    if sub == find(strcmp(settings.subjectDir,'s11')) && settings.cond == 5
        badTrials    = [1,badTrials]; %append virtual trial 1 to match EEG data
    end
    % Remove bad trials to align file with cleaned eeg data file
    A                = structfun(@(x) (removerows(x, 'ind', badTrials)), A, 'UniformOutput', false);
    T                = structfun(@(x) (removerows(x, 'ind', badTrials)), T, 'UniformOutput', false);

    % Identify all correct trials from that condition:
    corrTrial        = [];
    corrTrial        = find(A.points>=settings.pointsFilter & A.tempID==settings.cond & A.probeRT>0);

    % To get all incorrect trials
%     corrTrial1 = ones(1,dataLen)';
%     corrTrial1(corrTrial) = 0;
%     corrTrial = find(corrTrial1 == 1);
    
    % table containing variables for GLM
    glmTbl.(freq) = calc_glmTbl_timepoints('beh',settings,glmTbl,A,corrTrial,sub,freq,T);

    pressT = A.timing;
end