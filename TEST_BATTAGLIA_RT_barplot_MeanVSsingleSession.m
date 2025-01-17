%% function TEST_BATTAGLIA_RT_barplot_MeanVSsingleSession
clear; close all;
par.irng = 10;
rng(par.irng)
%% check if is on server
[~,nn]=system('hostname'); nn=strtrim(nn);
if strcmp(nn,'rk018610')
    isonserver=true;
else
    isonserver=false;
end
if ~isonserver
    test_dir    ='D:\SAPIENZA_Dataset\TEST_SAPIENZA';
else
    test_dir    ='~/SAPIENZA/SERVER/DCM';
end
%% Selezione dei file
directory_path = 'D:\SAPIENZA_Dataset\BATTAGLIA_Data';
all_files = dir(fullfile(directory_path, '*.mat'));
file_names = {all_files.name};
files_H = file_names(endsWith(file_names, 'H_Raw.mat'));
files_E = file_names(endsWith(file_names, 'E_Raw.mat'));

% fixed parameters
session_name = erase(files_H,'H_Raw.mat');
idir                            = 1;        % directions -> 1-8
S                               = filesep;
%% Step 0: arrange trials

% dir_path = strcat('D:\main_scriptDCM\','AllSessions','\rng',num2str(par.irng),filesep);


RT_session = struct();
lenTot = struct();
Bad_session = struct();
for idsess = 1:length(session_name)
    fprintf('Step 0: arrange trials\n');
    par.BattagliaArrangeTrials              = BattagliaArrangeTrialsParams();
    % par.BattagliaArrangeTrials.whichmodel   = 7;        % 7 Model for action session. 5 Model for all sessions
    par.BattagliaArrangeTrials.isdemo       = 1;        % getSelectionIndexes(1,1:3);
    par.BattagliaArrangeTrials.selS         = 1;
    par.BattagliaArrangeTrials.selK         = 1;
    par.BattagliaArrangeTrials.session_name = session_name{idsess};  % which session
    data_trials           = BattagliaArrangeTrials(par.BattagliaArrangeTrials);
    
    % Delete label 0 trials
    idx_empty = find(arrayfun(@(x) isempty(x.trialType), data_trials));
    
    if idx_empty ~= 0
        Bad_session(idsess).name = session_name(idsess);
        Bad_session(idsess).chamber = data_trials(idx_empty).Chamber;
        Bad_session(idsess).trialsId = data_trials(idx_empty).trialId;
        Bad_session(idsess).Trials = data_trials(idx_empty);
    else
        Bad_session(idsess).name = session_name(idsess);
        Bad_session(idsess).chamber = data_trials(idsess).Chamber;
        Bad_session(idsess).trialsId = [];
        Bad_session(idsess).Trials = [];
    end
    data_trials(idx_empty) = [];
    % add trialName
    [~,Labels] = getJointMonkeysLabels(1:24);
    for iTrial=1:length(data_trials)
        data_trials(iTrial).trialName        = Labels{data_trials(iTrial).trialType};
    end
    %% Step 1: Visual data
    data_trials = findSKdirection(data_trials);

    num_dir = 8; % direzioni movimento del task
    num_cond = 3; % numero condizioni 1-SoloS 2-SoloK 3-Joint S-K
    condition_name  = {'Act S-Obs K';'Obs S-Act K';'Act S-Act K'};
    color_name = {[0 0 0.5],[0 0.5 0],[1 0.6 0]};

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
    for cd=1:num_cond
        for ndir=1:num_dir
            Move_app = Move(cd).(strcat('dir',num2str(ndir)));
            RT_Sdir = NaN(lenM(cd,ndir),1);
            RT_Kdir = NaN(lenM(cd,ndir),1);
            for k = 1:lenM(cd,ndir)
                RT_Sdir(k,1) = Move_app(k).(strcat(RTfieldname,'_S'));
                RT_Kdir(k,1) = Move_app(k).(strcat(RTfieldname,'_K'));
            end
            RT_S(cd).(strcat('dir',num2str(ndir))) = RT_Sdir;
            RT_K(cd).(strcat('dir',num2str(ndir))) = RT_Kdir;
            RT_NaN(cd).(strcat('dir',num2str(ndir))).len = length(RT_Sdir);
            RT_NaN(cd).(strcat('dir',num2str(ndir))).S = [sum(isnan(RT_Sdir))];
            RT_NaN(cd).(strcat('dir',num2str(ndir))).K = [sum(isnan(RT_Kdir))];
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
    RT_session(idsess).cond = RT_bar;
    lenTot(idsess).len = lenM;
end
RT_total = struct();
nSession = numel(RT_session); 
for condIdx = 1:num_cond+1
    for dirIdx = 1:num_dir
        dirName = sprintf('dir%d', dirIdx);
        tempData = struct();
        for nSess = 1:nSession
            if isfield(RT_session(nSess).cond(condIdx), dirName)
                tempData(nSess).(dirName) = RT_session(nSess).cond(condIdx).(dirName);
            end
        end
        RT_total(condIdx).(dirName) = vertcat(tempData.(dirName));
    end
end
%% RT barplot
infoNaN = stringNaN(RT_NaN,num_cond,num_dir);
categories= {'dir1','dir2', 'dir3', 'dir4', 'dir5', 'dir6', 'dir7', 'dir8'};
RT_Mean = struct();
RT_Std = struct();
for icond = 1:num_cond+1
    for ndir = 1:ndir
        RT_calc = RT_total(icond).(strcat('dir',num2str(ndir)));
        RT_Mean(icond).(strcat('dir',num2str(ndir))) = 1000*mean(RT_calc);
        RT_Std(icond).(strcat('dir',num2str(ndir))) = 1000*std(RT_calc);
    end
end
RT_plotMean = (cell2mat(struct2cell(RT_Mean')))'; 
RT_plotStd = (cell2mat(struct2cell(RT_Std')))';
max_y = max(max(RT_plotMean+2*RT_plotStd));
%% BOXPLOT FIGURE

RT_Singlesession = struct();
lenTot = struct();
Bad_sessionSingle = struct();
session_name = 'SK009';
fprintf('Step 0: arrange trials\n');
par.BattagliaArrangeTrials              = BattagliaArrangeTrialsParams();
% par.BattagliaArrangeTrials.whichmodel   = 7;        % 7 Model for action session. 5 Model for all sessions
par.BattagliaArrangeTrials.isdemo       = 1;        % getSelectionIndexes(1,1:3);
par.BattagliaArrangeTrials.selS         = 1;
par.BattagliaArrangeTrials.selK         = 1;
par.BattagliaArrangeTrials.session_name = session_name;  % which session
data_trialsSingle           = BattagliaArrangeTrials(par.BattagliaArrangeTrials);
    
% Delete label 0 trials
idx_empty = find(arrayfun(@(x) isempty(x.trialType), data_trialsSingle));

if idx_empty ~= 0
    Bad_sessionSingle.name = session_name;
    Bad_sessionSingle.chamber = data_trialsSingle(idx_empty).Chamber;
    Bad_sessionSingle.trialsId = data_trialsSingle(idx_empty).trialId;
    Bad_sessionSingle.Trials = data_trialsSingle(idx_empty);
else
    Bad_sessionSingle.name = session_name;
    Bad_sessionSingle.chamber = data_trialsSingle.Chamber;
    Bad_sessionSingle.trialsId = [];
    Bad_sessionSingle.Trials = [];
end
data_trialsSingle(idx_empty) = [];
% add trialName
[~,LabelsSingle] = getJointMonkeysLabels(1:24);
for iTrial=1:length(data_trialsSingle)
    data_trialsSingle(iTrial).trialName        = LabelsSingle{data_trialsSingle(iTrial).trialType};
end
data_trialsSingle = findSKdirection(data_trialsSingle);

par.RTsingleEXTRACT.num_cond = num_cond;
par.RTsingleEXTRACT.num_dir = num_dir;
[RT_Ssingle,RT_KSingle,RTSingle_plotMean,RTSingle_plotStd,max_ySingle,RTSingle_name] = RTsingleEXTRACT(data_trialsSingle,par.RTsingleEXTRACT);

dir_path = strcat('D:\main_scriptDCM\',session_name,'\rng',num2str(par.irng),'\');

%% Barplot Figure
par.RTbarplotAllvsSingleDir.color_name = color_name;
par.RTbarplotAllvsSingleDir.num_dir = num_dir;
par.RTbarplotAllvsSingleDir.num_cond = num_cond;
par.RTbarplotAllvsSingleDir.max_y = max_y;
par.RTbarplotAllvsSingleDir.idcond = 1;
par.RTbarplotAllvsSingleDir.session_name = session_name;
par.RTbarplotAllvsSingleDir.chamber = data_trialsSingle(1).Chamber;
par.RTbarplotAllvsSingleDir.dir_path = dir_path;

RTbarplotAllvsSingleDir(RT_plotMean,RT_plotStd,RT_Ssingle,RT_KSingle,RTSingle_name,par.RTbarplotAllvsSingleDir);

close all;
