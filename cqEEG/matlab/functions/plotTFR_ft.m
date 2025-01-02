function plotTFR_ft(what,settings,TFR,freq,method,plottype,sub)

switch(what)
    case 'subTFR'
        switch(plottype)
            case 'singleplot'
                % plot TFR at C3
                cfg = [];
                cfg.xlim         = [-0.5 4.5];
                cfg.ylim         = [6 40];
                cfg.zlim         = [-30 30];
                cfg.channel      = 'D19'; % C3
                cfg.layout       = 'biosemi128.lay';

                fig = figure();
                ft_singleplotTFR(cfg, TFR);
                hold on;
                plot([0, 0], [6, 40], 'k', 'LineWidth', 1); % Add vertical line at 0
                plot([1.5, 1.5], [6, 40], '--k', 'LineWidth', 1); % Add dotted vertical line at 1.5
                hold off;
%                 title(sprintf('TFR %s %s %s %s %s',method,freq,cfg.channel,getCondStr(settings.cond),settings.subjectDir{sub}));

                %                 figure();
                %                 histogram(TFR.powspctrm(111,:,30:45));
                %                 hold on;
                %                 meanVal = nanmean(TFR.powspctrm(111,:,30:45),'all'); % baseline mean
                %                 line([meanVal, meanVal], ylim, 'Color', 'r', 'LineStyle', '-');
                %                 hold off;


                % 3D plot--------------------------------------------------
%                 figure();
%                 chanData = squeeze(TFR.powspctrm(:,111,:,:));
%                 timeAxis = TFR.time;
%                 freqAxis = TFR.freq;
%                 [time, freq] = meshgrid(timeAxis, freqAxis);
%                 surf(time, freq, squeeze(mean(chanData, 1)), 'EdgeColor', 'none');
%                 title('3D TFR Heatmap for Channel 111');
%                 xlabel('Time (s)');
%                 ylabel('Frequency (Hz)');
%                 zlabel('Power (% change)');
% -------------------------------------------------------------------------

            case 'topoplot'
                for period = settings.period
                    period=period{:};
                    cfg              = [];
                    cfg.xlim         = settings.(getCondStr(settings.cond)).(period); % time interval: preparation period
                    %                     cfg.zlim         = [-Inf 0.1]; % adjust scale - interested in power changes <0 -> why???
                    cfg.colorbar     = 'yes';
                    cfg.marker       = 'on';
                    cfg.layout       = 'biosemi128.lay';
                    cfg.comment      = 'no';   % Turn off any default commentary in the plot

                    fig = figure();
                    ft_topoplotTFR(cfg, TFR);
%                     title(sprintf('TFR %s %s %s %s %s',method,freq,getCondStr(settings.cond),settings.subjectDir{sub},period));

                    if settings.savePlt
                        condStr = getCondStr(settings.cond);
                        filepath = fullfile(settings.saveDirTFRvFigs, 'sub/', settings.subjectDir{sub});
                        checkDir(filepath)
                        saveas(fig,fullfile(filepath, sprintf('%sTFR%s_%s_%s_%s_s%02i',plottype,method,freq,condStr,period,settings.subjectnumber(sub))), 'fig');
                    end
                end
        end

    case 'gaTFR'
        switch(plottype)
            case 'singleplot'
                cfg              = [];
                cfg.xlim         = [-0.5 4.5];
                cfg.ylim         = [6 40];
                cfg.zlim         = [-30 30];
                cfg.channel      = 'all'; %C3-contralateral M1- Consider neighbours? - the default takes avg power across all channels
                cfg.layout       = 'biosemi128.lay'; % Layout for the Biosemi high density EEG system @ Bangor
                cfg.comment = 'no';   % Turn off any default commentary in the plot

                fig = figure();
                ft_singleplotTFR(cfg, TFR);
                hold on;
                plot([0, 0], [6, 40], 'k', 'LineWidth', 1); % Add vertical line at 0
                plot([1.5, 1.5], [6, 40], '--k', 'LineWidth', 1); % Add dotted vertical line at 1.5
                hold off;
                title('');
%                 title(sprintf('TFR %s %s %s',method,cfg.channel,getCondStr(settings.cond)));

            case 'topoplot'
                for period = settings.period
                    period=period{:};
                    cfg              = [];
                    cfg.xlim         = settings.(getCondStr(settings.cond)).(period); % time interval: preparation period
                    %                     cfg.zlim         = [-.2 0]; % adjust scale - interested in power changes <0 -> why???
                    cfg.colorbar     = 'yes';
                    cfg.marker       = 'on';
                    cfg.layout       = 'biosemi128.lay';
                    cfg.comment      = 'no';   % Turn off any default commentary in the plot

                    fig = figure();
                    ft_topoplotTFR(cfg, TFR);
                    title('');
%                     title(sprintf('TFR %s %s %s %s',method,freq,getCondStr(settings.cond),period));

                    if settings.savePlt
                        condStr = getCondStr(settings.cond);
                        filepath = fullfile(settings.saveDirTFRvFigs, 'group/');
                        checkDir(filepath)
                        saveas(fig,fullfile(filepath, sprintf('%sTFR%s_%s_%s_%s',plottype,method,freq,condStr,period)), 'fig');
                    end
                end

            case 'traceplot'
                fig = figure();
                ax = axes('Parent',fig);
                hold(ax,'on')
                % squeeze(nanmean(TFR.powspctrm(settings.chan,:,:),2))
                meanBeta = nanmean(TFR.powspctrm,2);
                meanBeta = squeeze(meanBeta);
                timeAx =  TFR.time;
                smoothBins=2;
                smoothType=1;
                smoothEnds=1;
                smoothYerr  = fastsmooth(nanstd(meanBeta(:,:),[],1),smoothBins,smoothType,smoothEnds);
                for i = 1:size (meanBeta,2)
                    meanBetaSmoothed(:,i) = fastsmooth(meanBeta(:,i),smoothBins,smoothType,smoothEnds);
                end

                % shadedErrorBar(timeAx,mean(meanBetaSmoothed(:,:),1),smoothYerr./sqrt(size(meanBeta,1)),'lineprops',{'color',[0 0 0]}); % sem
                shadedErrorBar(timeAx,mean(meanBetaSmoothed(:,:),1),std(meanBeta),'lineprops','b'); % std

                axis square
                set(gca,'FontSize',10)

                axis([-.5 4.5 -50 50])
                plot([0, 0], [-50, 50], 'k', 'LineWidth', 1); % Add vertical line at 0
                plot([1.5, 1.5], [-50, 50], '--k', 'LineWidth', 1); % Add dotted vertical line at 1.5
                %             legend('Beta band')
                xlabel('Time (s)')
                ylabel('Relative change (%)') % ylabel('Power (\mu V^2)')
                title(sprintf('gaTFR traceplot %s %s',getCondStr(settings.cond),freq));
                hold off

            case 'topoplotMarkers'

                load('/rds/projects/k/kornyshk-kornyshevalab/martin/cqEEG/data/eeg43/tfr/voltage/data/ana/erdChans.mat')
                filteredTable = maxERDchannels((maxERDchannels.Condition == settings.cond) & strcmp(maxERDchannels.Frequency, freq), :);
                uniqueLabels = unique(filteredTable.Label);
                indicesArray = [];
                for i = 1:length(uniqueLabels)
                    label = uniqueLabels{i};
                    labelIdx = find(strcmp(TFR.label, label));
                    indicesArray = [indicesArray; labelIdx]; % Append indices to the array
                end

                % Create a dummy TFR structure with sensor labels and no data
                dummyTFR = [];
                dummyTFR.label = TFR.label;
                dummyTFR.dimord = 'chan_freq_time';
                dummyTFR.freq = 1;
                dummyTFR.time = 1;
                dummyTFR.powspctrm = zeros(length(dummyTFR.label), 1, 1);  % No power data, just zeros
                
                % Topoplot config
                cfg = [];
                cfg.parameter = 'powspctrm';
                cfg.xlim = [1 1];
                cfg.marker = 'on';
                cfg.layout = 'biosemi128.lay';
                cfg.highlight = 'on';
                cfg.highlightchannel = {[indicesArray]}; % C3
%                 cfg.highlightchannel = {[110,113,120,122]}; % Cluster 4
%                 cfg.highlightchannel = {[96,104,106,108,110,113,115,120,122]}; % Cluster 9
                cfg.highlightsymbol = '.';
                cfg.highlightcolor = {'r'};
                cfg.highlightsize = [40];
                cfg.colorbar = 'no'; % No colorbar needed
                cfg.style = 'blank';
                cfg.comment = 'no';   % Turn off any default commentary in the plot
                
                % Create the figure
                fig = figure();
                ft_topoplotTFR(cfg, dummyTFR);
                title(sprintf('Sensor Position Plot %s %s',getCondStr(settings.cond),freq));
                return
        end
end

if settings.savePlt
    condStr = getCondStr(settings.cond);
    switch(what)
        case 'subTFR'
            filepath = fullfile(settings.saveDirTFRvFigs, 'sub/', settings.subjectDir{sub});
            checkDir(filepath)
            switch(plottype)
                case 'singleplot'
                    saveas(fig,fullfile(filepath, sprintf('%sTFR%s_%s_%s_s%02i',plottype,method,freq,condStr,settings.subjectnumber(sub))), 'fig');
                case 'topoplot'
                    % easiest to save it above already
            end
        case 'gaTFR'
            filepath = fullfile(settings.saveDirTFRvFigs, 'group/');
            checkDir(filepath)
            switch(plottype)
                case 'singleplot'
                    saveas(fig,fullfile(filepath, sprintf('gaTFR%s_%s_%s',method,freq,condStr)), 'fig');
                case 'traceplot'
                    saveas(fig,fullfile(filepath, sprintf('trace_gaTFR%s_%s_%s',method,freq,condStr)), 'fig');
            end
    end
end
end