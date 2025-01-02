function [velocity] = calc_TFRvel(TFR)

    % derivative of TFR
    derivativeTFR = diff(TFR.powspctrm, 1, 4);

    % calculate velocity of signal at go cue
    velocity.goCue = squeeze(nanmean(derivativeTFR(:,111,:,29),3));

    % calculate mean velocity over 3 time points ending with go cue
    velocity.hiLo = squeeze(nanmean(nanmean(derivativeTFR(:,111,:,27:29),3),4));

    % plot TFR and its derivative (=velocity)
%     pltTrial = 13;
%     figure();
%     pltTFR = plot(squeeze(nanmean(TFR.powspctrm(pltTrial,111,:,:),3)), 'DisplayName', 'TFR');
%     hold on;
%     pltVel = plot(squeeze(nanmean(derivativeTFR(pltTrial,111,:,:),3)), 'DisplayName', 'velocity');
%     legend([pltTFR, pltVel], 'Location', 'northeast');
%     hold off;

end