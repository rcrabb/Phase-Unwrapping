function avebyfield(s)

sfields = fieldnames(s);
slen = length(s);

for i = 1:length(sfields)
    fname = sfields{i};
    if strcmp(fname,'pctcorrect') || strcmp(fname,'States') || strcmp(fname,'m') || strcmp(fname,'endstate') 
        continue
    end
    if isempty(getfield(s,{1},fname))
        continue
    end
    disp(fname);
    % Wow can't seem to get all the items from a field at once
    for j = 1:slen
        tmp(j) = getfield(s,{j},fname);
    end
    uniques = unique(tmp);
    clear tmp;
    if length(uniques>1)
    for u = uniques
        cnt = 0;
        for j = 1:slen
            if getfield(s,{j},fname) == u
                cnt = cnt+1;
                tmp(cnt) = getfield(s,{j},'pctcorrect');
            end
        end
        disp(['Bin ' num2str(u) ': | mean:' num2str(mean(tmp)) ' | median:' num2str(median(tmp)) ' | max:' num2str(max(tmp))]);
    end
    end
    
    
end