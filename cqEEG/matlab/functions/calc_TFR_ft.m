function [TFR] = calc_TFR_ft(what,settings,data,freq,corrTrial,method)

% switch(what)
%     case 'subTFR'
        switch(method)
            case 'Morl'
                % calculate TFR: Morlet wavelet
                cfg              = [];
                cfg.output       = 'pow';
                cfg.channel      = 'all';
                cfg.trials       = corrTrial; % ft_freqanalysis will read in all trials from level data but will analyse only trlNr
                cfg.method       = 'wavelet';
                cfg.width        = 7; 
                cfg.toi          = -1.3:.1:5.3; % in sliding windows of 0.1 s     
                cfg.foi          = settings.(freq).foi;
                cfg.keeptrials   = 'yes';
                cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.1; % length of time window = 0.1 sec
                
                TFR              = ft_freqanalysis(cfg, data); % TFR.powspctrm contains the power for chanxfreqxtime

            case 'Hann'
                % TFR: Hanning taper, fixed window length
                cfg              = [];
                cfg.output       = 'pow';
                cfg.channel      = 'all';
                cfg.trials       = corrTrial;
                cfg.method       = 'mtmconvol';
                cfg.taper        = 'hanning';
                cfg.foi          = settings.(freq).foi;
                cfg.keeptrials   = 'yes';
                cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.1; % length of time window = 0.1 sec
                cfg.toi          = data.time{1,1}(:,[1:100:find(data.time{1,1}(:,:)>=4.5+0.8,1,'first')]); % in sliding windows of 0.1 s     
                
                TFR              = ft_freqanalysis(cfg, data); % TFRhann.powspctrm contains the power for chnxfreqxtime

        end
% end
end






