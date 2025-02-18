
function [VEL_resultACT,VEL_resultJOINT,VEL_countJOINT,VEL_percJOINT]=countVELoccurrance(VEL_total,par)

num_dir = length(fields(VEL_total));

for dirIdx=1:num_dir
    dirName = sprintf('dir%d', dirIdx);
    VEL_ActS = VEL_total(1).(dirName);
    VEL_ActK = VEL_total(2).(dirName);
    VEL_S = VEL_total(3).(dirName);
    VEL_K = VEL_total(4).(dirName);
    if mean(VEL_ActS) > mean(VEL_ActK)
        VEL_resultACT(1).(dirName) = 'S';
    else
        VEL_resultACT(1).(dirName) = 'K';
    end
    count_S = 0;
    count_K = 0;
    for i=1:length(VEL_S)
        if VEL_S(i) > VEL_K(i)
            count_S = count_S+1;
        else
            count_K = count_K +1;
        end
    end
    VEL_countJOINT(1).(dirName) = [count_S;count_K];
    VEL_percJOINT(1).(dirName) = [100*count_S/length(VEL_S);100*count_K/length(VEL_K)];
end


for dirIdx=1:num_dir
    dirName = sprintf('dir%d', dirIdx);
    S = VEL_percJOINT(1).(dirName)(1);
    K = VEL_percJOINT(1).(dirName)(2);
    if S>K
        if S>par.perc
            VEL_resultJOINT(1).(dirName) = 'S';
        else
            VEL_resultJOINT(1).(dirName) = '50%';
        end
    else
        if K>par.perc
            VEL_resultJOINT(1).(dirName) = 'K';
        else
            VEL_resultJOINT(1).(dirName) = '50%';
        end
    end
end

