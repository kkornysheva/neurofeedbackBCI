function [TFR] = rem_powOutliers(TFR)

    % outlier removal: 3 Standard Deviations
    meanVal = nanmean(TFR.powspctrm(:));
    stdVal = nanstd(TFR.powspctrm(:));
    
    outlierIdx = abs(TFR.powspctrm - meanVal) > (3 * stdVal);
    TFR.powspctrm(outlierIdx) = NaN;
    
%     numNotNaN = sum(~isnan(TFR.powspctrm(:)));
end