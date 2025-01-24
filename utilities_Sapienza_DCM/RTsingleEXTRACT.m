
function [RT_S,RT_K,RT_plotMean,RT_plotStd,max_y,RT_name] = RTsingleEXTRACT(data_trials,par)

num_cond = par.num_cond;
num_dir = par.num_dir;

Condition = cell(num_cond,num_dir);

for cd=1:num_cond
    for ndir=1:num_dir
        Condition{cd,ndir} = find([data_trials.Condition]==cd & [data_trials.Direction]==ndir);
    end
end

Move = struct();
lenM = NaN(num_cond,num_dir);
for cd=1:num_cond
    for ndir=1:num_dir
        Move(cd).(strcat('dir',num2str(ndir))) = data_trials(Condition{cd,ndir});
        lenM(cd,ndir) = length(data_trials(Condition{cd,ndir}));
    end
end
%% Reaction Time
RTfieldname = 'RT';
RT_S = struct();
RT_K = struct();
% RT_Smean = struct();
% RT_Kmean = struct();
% RT_Sstd = struct();
% RT_Kstd = struct();
RT_NaN = struct();
RT_name = struct();
for cd=1:num_cond
    for ndir=1:num_dir
        Move_app = Move(cd).(strcat('dir',num2str(ndir)));
        RT_Sdir = NaN(lenM(cd,ndir),1);
        RT_Kdir = NaN(lenM(cd,ndir),1);
        RT_trialID = NaN(lenM(cd,ndir),1);
        for k = 1:lenM(cd,ndir)
            RT_Sdir(k,1) = Move_app(k).(strcat(RTfieldname,'_S'));
            RT_Kdir(k,1) = Move_app(k).(strcat(RTfieldname,'_K'));
            RT_trialID(k,1)   = Move_app(k).trialId;
        end
        RT_S(cd).(strcat('dir',num2str(ndir))) = RT_Sdir;
        RT_K(cd).(strcat('dir',num2str(ndir))) = RT_Kdir;
        RT_NaN(cd).(strcat('dir',num2str(ndir))).len = length(RT_Sdir);
        RT_NaN(cd).(strcat('dir',num2str(ndir))).S = [sum(isnan(RT_Sdir))];
        RT_NaN(cd).(strcat('dir',num2str(ndir))).K = [sum(isnan(RT_Kdir))];
        RT_name(cd).(strcat('dir',num2str(ndir))).D = ndir;
        RT_name(cd).(strcat('dir',num2str(ndir))).T = RT_trialID;
        %     RT_Smean(cd).(strcat('dir',num2str(ndir))) = mean(RT_Sdir,'omitmissing');
        %     RT_Kmean(cd).(strcat('dir',num2str(ndir))) = mean(RT_Kdir,'omitmissing');
        %     RT_Sstd(cd).(strcat('dir',num2str(ndir))) = std(RT_Sdir,'omitmissing');
        %     RT_Kstd(cd).(strcat('dir',num2str(ndir))) = std(RT_Kdir,'omitmissing');
    end
end

%% Braplot per DIREZIONE
RT_bar(1) = RT_S(1);
RT_bar(2) = RT_K(2);
RT_bar(3) = RT_S(3);
RT_bar(4) = RT_K(3);

%% RT barplot
RT_Mean = struct();
RT_Std = struct();
for icond = 1:num_cond+1
    for ndir = 1:ndir
        RT_calc = RT_bar(icond).(strcat('dir',num2str(ndir)));
        RT_Mean(icond).(strcat('dir',num2str(ndir))) = 1000*mean(RT_calc);
        RT_Std(icond).(strcat('dir',num2str(ndir))) = 1000*std(RT_calc);
    end
end
RT_plotMean = (cell2mat(struct2cell(RT_Mean')))';
RT_plotStd = (cell2mat(struct2cell(RT_Std')))';
max_y = max(max(RT_plotMean+2*RT_plotStd));