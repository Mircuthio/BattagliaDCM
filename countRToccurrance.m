
function [RT_resultACT,RT_resultJOINT,RT_countJOINT,RT_percJOINT]=countRToccurrance(RT_total,par)

num_dir = length(fields(RT_total));
for dirIdx=1:num_dir
    dirName = sprintf('dir%d', dirIdx);
    RT_ActS =RT_total(1).(dirName);
    RT_ActK =RT_total(2).(dirName);
    RT_S = RT_total(3).(dirName);
    RT_K = RT_total(4).(dirName);
    if mean(RT_ActS) < mean(RT_ActK)
       RT_resultACT(1).(dirName) = 'S';
    else
       RT_resultACT(1).(dirName) = 'K';
    end
    count_S = 0;
    count_K = 0;
    for i=1:length(RT_S)
        if RT_S(i)<RT_K(i)
            count_S = count_S+1;
        else
            count_K = count_K +1;
        end
    end
    RT_countJOINT(1).(dirName) = [count_S;count_K];
    RT_percJOINT(1).(dirName) = [100*count_S/length(RT_S);100*count_K/length(RT_K)];
end


for dirIdx=1:num_dir
    dirName = sprintf('dir%d', dirIdx);
    S = RT_percJOINT(1).(dirName)(1);
    K = RT_percJOINT(1).(dirName)(2);
    if S>K
        if S>par.perc
            RT_resultJOINT(1).(dirName) = 'S';
        else
            RT_resultJOINT(1).(dirName) = '50%';
        end
    else
        if K>par.perc
            RT_resultJOINT(1).(dirName) = 'K';
        else
            RT_resultJOINT(1).(dirName) = '50%';
        end
    end
end

