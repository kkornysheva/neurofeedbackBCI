% GLMs and LMEs

cd /rds/projects/k/kornyshk-kornyshevalab/martin/cqEEG/

close all; clear; clc;

addpath(genpath('./functions'));

set.dependentVar = 'RT'; % 'RT' or 'points' or 'timingPrec'
set.freq = {'theta','mu','lowBeta','highBeta','beta','all','maxERD'};
% set.freq = {'theta'};
set.condID = [1,2,5];

set.powPredictors = {'power','pow0_300','pow400_700','pow800_1100','pow1200_1500','prodStart','velGo','velHiLo'};
% set.powPredictors = {'power0_300','pow300_600','pow600_900','pow900_1200','pow1200_1500','prodStart','velGo','velHiLo'};

load('./data/eeg63/data/glm/glmTbl.mat')

% subjectMaxERD = unique(glmTbl.maxERD(:, {'subject', 'maxERDfreq'}), 'rows');
% mean(subjectMaxERD.maxERDfreq);
% std(subjectMaxERD.maxERDfreq);

% figure();histogram(glmTbl.mu.pow0_300); hold on; line([median(glmTbl.mu.pow0_300), median(glmTbl.mu.pow0_300)], ylim, 'Color', 'r', 'LineWidth', 2);hold off;

%% LME: Trial-by-trial

resultsStruct.Intercept = struct('Frequency', {}, 'Power', {}, 'CoefName', {}, 'Estimate', {}, 'SE', {}, 'PValue', {}, 'CI', {});
resultsStruct.Predictor1 = struct('Frequency', {}, 'Power', {}, 'CoefName', {}, 'Estimate', {}, 'SE', {}, 'PValue', {}, 'CI', {});
resultsStruct.Predictor2 = struct('Frequency', {}, 'Power', {}, 'CoefName', {}, 'Estimate', {}, 'SE', {}, 'PValue', {}, 'CI', {});

% Loop through frequencies and power predictors
for freq = set.freq
    freq = freq{:};
    
    for pow = set.powPredictors
        pow = pow{:};

        % Remove outliers above 100 ---------------------------------------
%         glmTbl.(freq).(pow)(glmTbl.(freq).(pow) > 100) = NaN;
%         sglmTbl.(freq).RT(glmTbl.(freq).RT < 150) = NaN;
        % -----------------------------------------------------------------

        % Specify formula
        formula = sprintf('%s ~ 1 + %s  + (1 + %s |subject)', set.dependentVar, pow, pow);
%         formula = sprintf('%s ~ 1 + %s  + condition + (1 + %s |subject) + (1 + condition |subject)', set.dependentVar, pow, pow);

        % Fit the linear mixed-effects model
        lmeModel = fitlme(glmTbl.(freq), formula);

        % Extract coefficients, p-values, and confidence intervals
        coefTable = table(lmeModel.Coefficients.Name, lmeModel.Coefficients.Estimate, lmeModel.Coefficients.SE, lmeModel.Coefficients.pValue, 'VariableNames', {'Name', 'Estimate', 'SE', 'PValue'});
        
        % Extract confidence intervals
        CI = coefCI(lmeModel);
        coefTable.CI = string([CI(:, 1) CI(:, 2)]);

        % Add results to the struct
        resultsStruct.Intercept(end+1) = struct('Frequency', freq, 'Power', pow, 'CoefName', coefTable.Name(1), 'Estimate', coefTable.Estimate(1), 'SE', coefTable.SE(1), 'PValue', coefTable.PValue(1), 'CI', CI(1, :));
        resultsStruct.Predictor1(end+1) = struct('Frequency', freq, 'Power', pow, 'CoefName', coefTable.Name(2), 'Estimate', coefTable.Estimate(2), 'SE', coefTable.SE(2), 'PValue', coefTable.PValue(2), 'CI', CI(2, :));
%         resultsStruct.Predictor2(end+1) = struct('Frequency', freq, 'Power', pow, 'CoefName', coefTable.Name(3), 'Estimate', coefTable.Estimate(3), 'SE', coefTable.SE(3), 'PValue', coefTable.PValue(3), 'CI', CI(3, :));
    end
end

%% GLM per subject incl. visualisation

for freq = set.freq
    freq=freq{:};
    
    % Create a figure for each frequency
    figure;
    hold on;

    % Plot a scatter plot for each subject
    subjects = unique(glmTbl.(freq).subject);
    for subjIdx = 1:length(subjects)
        subject = subjects(subjIdx);
       
        subjectData = glmTbl.(freq)(glmTbl.(freq).subject == subject, :);

        % Remove outliers above 100 ---------------------------------------
        for j = 1:length(set.powPredictors)
            currPred = set.powPredictors{j};
            currCol = subjectData.(currPred);
            currCol(currCol > 100) = NaN;
            subjectData.(currPred) = currCol;
        end
        % -----------------------------------------------------------------
        % Remove RTs below 150 ms
%         subjectData.RT(subjectData.RT < 150) = NaN;

        subplot(5,4,subject);
        hold on;
        switch(set.dependentVar)
            case 'RT'
                scatter(subjectData.(set.powPredictors{2}), subjectData.RT, 'DisplayName', ['Subject ' num2str(subject)]);%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'points'
                scatter(subjectData.power, subjectData.points, 'DisplayName', ['Subject ' num2str(subject)]); 
            case 'timingPrec'
                scatter(subjectData.power, subjectData.timingPrec, 'DisplayName', ['Subject ' num2str(subject)]); 
        end

        % Specify formula
        formula = sprintf('%s ~ 1 + %s', set.dependentVar, set.powPredictors{2}); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Fit a linear regression model for each subject using fitglm
        model = fitglm(subjectData, formula);
        
        % Predicted values using the fitted model
        yFit = predict(model, subjectData);
        
        % Plot the regression line for each subject
        
        plot(subjectData.(set.powPredictors{2}), yFit, 'LineWidth', 2, 'DisplayName', ['Regression Line - Subject ' num2str(subject)]); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Add histogram to the subplot
        histogram(subjectData.(set.powPredictors{2}), 30, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1.5); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        title(['subject', num2str(subject)]);
        legend('Location', 'northeast');
        legend(['Estimate: ', num2str(model.Coefficients{'pow0_300','Estimate'}),'p-val: ', num2str(model.Coefficients{'pow0_300','pValue'})]);
        
        pause(0.01);
    end

    title(freq);
    xlabel('Power');
    ylabel('RT');
    
    hold off;
%      figure();histogram(subjectData.(set.powPredictors{3}),30);
end

%% GLM per subject for plotting

for freq = set.freq
    freq=freq{:};
    
    % Create a figure for each frequency
    figHandle = figure();
    hold on;

    % Plot a scatter plot for each subject
    subjects = unique(glmTbl.(freq).subject);
    for subjIdx = 3
        subject = subjects(subjIdx);
       
        subjectData = glmTbl.(freq)(glmTbl.(freq).subject == subject, :);

        % Remove outliers above 100 ---------------------------------------
        for j = 1:length(set.powPredictors)
            currPred = set.powPredictors{j};
            currCol = subjectData.(currPred);
            currCol(currCol > 100) = NaN;
            subjectData.(currPred) = currCol;
        end
        % -----------------------------------------------------------------
        % Remove RTs below 150 ms
%         subjectData.RT(subjectData.RT < 150) = NaN;

%         subplot(5,4,subject);
        hold on;
        switch(set.dependentVar)
            case 'RT'
                scatter(subjectData.(set.powPredictors{2}), subjectData.RT, 'DisplayName', ['Subject ' num2str(subject)]);%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'points'
                scatter(subjectData.power, subjectData.points, 'DisplayName', ['Subject ' num2str(subject)]); 
            case 'timingPrec'
                scatter(subjectData.power, subjectData.timingPrec, 'DisplayName', ['Subject ' num2str(subject)]); 
        end

        % Specify formula
        formula = sprintf('%s ~ 1 + %s', set.dependentVar, set.powPredictors{2}); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Fit a linear regression model for each subject using fitglm
        model = fitglm(subjectData, formula);
        
        % Predicted values using the fitted model
        yFit = predict(model, subjectData);
        
        % Plot the regression line for each subject
        
        plot(subjectData.(set.powPredictors{2}), yFit, 'LineWidth', 2, 'DisplayName', ['Regression Line - Subject ' num2str(subject)]); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Add histogram to the subplot
%         histogram(subjectData.(set.powPredictors{3}), 30, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1.5); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    hold off;
    exportgraphics(figHandle,'./figures/GLM_sub3.eps','ContentType','vector','BackgroundColor','none');

%      figure();histogram(subjectData.(set.powPredictors{3}),30);
end

%% LME: Trial-by-trial: multiple predictors

% Initialize results struct
resultsStruct.Intercept = struct('Frequency', {}, 'Power', {}, 'CoefName', {}, 'Estimate', {}, 'SE', {}, 'PValue', {}, 'CI', {});
resultsStruct.Predictor = table();


% Loop through frequencies
for freq = set.freq
    freq = freq{:};
    
    % Create a formula with all power predictors
    formula = sprintf('%s ~ %s + %s + %s + %s + (1|subject)', set.dependentVar, set.powPredictors{2:end});

    % Fit the linear mixed-effects model
    lmeModel = fitlme(glmTbl.(freq), formula);

    % Extract coefficients, p-values, and confidence intervals
    coefTable = table(lmeModel.Coefficients.Name, lmeModel.Coefficients.Estimate, lmeModel.Coefficients.SE, lmeModel.Coefficients.pValue, 'VariableNames', {'Name', 'Estimate', 'SE', 'PValue'});
    
    % Extract confidence intervals
    CI = coefCI(lmeModel);
    coefTable.CI = string([CI(:, 1) CI(:, 2)]);

    % Add results to the struct
    interceptResults = struct('Frequency', freq, 'Power', 'Intercept', 'CoefName', coefTable.Name(1), 'Estimate', coefTable.Estimate(1), 'SE', coefTable.SE(1), 'PValue', coefTable.PValue(1), 'CI', CI(1, :));

    % Loop through predictors and fill predictorResults struct
    predictorResults = table();
    for i = 2:numel(powPredictors) + 1
        predictorName = powPredictors{i - 1};
        predictorIdx = find(contains(coefTable.Name, predictorName));
    
        predictorResults = [predictorResults; table({freq}, {predictorName}, coefTable.Name(predictorIdx), coefTable.Estimate(predictorIdx), coefTable.SE(predictorIdx), coefTable.PValue(predictorIdx), CI(predictorIdx, 1), CI(predictorIdx, 2), 'VariableNames', {'Frequency', 'Power', 'CoefName', 'Estimate', 'SE', 'PValue', 'CI_Lower', 'CI_Upper'})];
    end

    resultsStruct.Intercept = [resultsStruct.Intercept; interceptResults];
    resultsStruct.Predictor = [resultsStruct.Predictor; predictorResults];
end

% Display the final results
% disp('Final Results:');
% disp(resultsStruct);

%% Plot histogram of RTs for conditions to check for condition differences in RT

uniqueCond = unique(glmTbl.lowBeta.condition);
ploti = 1;

figure;

% Loop through each condition and plot histogram
for cond = set.condID(1:end)
    disp(cond)
    conditionData = glmTbl.lowBeta(glmTbl.lowBeta.condition == cond, :);

    % Plot histogram
    subplot(1, length(uniqueCond), ploti);
    histogram(conditionData.RT);
    title(sprintf('Condition %s',getCondStr(cond)));

    % Add vertical line for median RT
    hold on;
    medianRT = nanmedian(conditionData.RT);
    xline(medianRT,'--r', 'LineWidth', 2, 'Label', ['Median RT: ' num2str(medianRT)], 'DisplayName', 'Median RT');
    hold off;
    ploti = ploti + 1;
end

%% Look at correct vs incorrect trials

set.dependentVar = 'correct';
corr = load('./data/eeg24/tfr/voltage/data/glm/glmTbl.mat')
incorr = load('./data/eeg30/tfr/voltage/data/glm/glmTbl.mat')
for freq = set.freq
    freq = freq{:};

    currentTbl = corr.glmTbl.(freq);
    currentTbl.correct = ones(height(currentTbl),1);
    corr.glmTbl.(freq) = currentTbl;

    currentTbl = incorr.glmTbl.(freq);
    currentTbl.correct = zeros(height(currentTbl),1);
    incorr.glmTbl.(freq) = currentTbl;

    glmTbl.(freq) = [corr.glmTbl.(freq); incorr.glmTbl.(freq)];
end
