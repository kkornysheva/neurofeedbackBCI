function checkDir(filepath)
    if ~isdir(filepath)
        mkdir(filepath);
    end
end