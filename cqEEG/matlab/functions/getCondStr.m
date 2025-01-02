function [condStr] = getCondStr(cond)
    if cond == 1
        condStr = "slow";
    elseif cond ==2
        condStr = "fast";
    elseif cond == 5
        condStr = "thumb";
    end   
end

