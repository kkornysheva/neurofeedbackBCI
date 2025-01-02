function saveTFR(what,settings,sub,method,TFR,freq)
    if ~settings.saveDat
        return
    end

    condStr = getCondStr(settings.cond);
    switch(what)
        case 'subTFR'
            filepath = fullfile(settings.saveDirTFRvData, 'sub/', settings.subjectDir{sub});
            checkDir(filepath)
            save(fullfile(filepath, sprintf('TFR%s_%s_%s_%s.mat',method,freq,condStr,settings.subjectDir{sub})),'TFR');
        case 'gaTFR'
            filepath = fullfile(settings.saveDirTFRvData, 'group/');
            checkDir(filepath)
            save(fullfile(filepath, sprintf('gaTFR%s_%s_%s.mat',method,freq,condStr)),'TFR');
    end
end