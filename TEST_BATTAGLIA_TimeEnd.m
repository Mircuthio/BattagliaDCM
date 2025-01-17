%% function TEST_BATTAGLIA_TimeEnd.m
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
    filesep = '\';
else
    test_dir    = '~/SAPIENZA/SERVER/DCM';
    filesep = '/';

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

dir_path = strcat('D:\main_scriptDCM\','AllSessions','\rng',num2str(par.irng),filesep);

ET_subj = struct();
ET_subRis =struct();
ET_subMax =struct();
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
    data_trials                             = BattagliaArrangeTrials(par.BattagliaArrangeTrials);
    
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

    %% Cursor Movement (x,y)
    Jfieldname = 'JXYEXY';
    Jtime = 'timeET';
    tStart = -0.2;
    ET_session = struct();
    ET_rispSess = struct();
    ET_maxSess = struct();
    for cd=1:num_cond

        for ndir=1:num_dir
            M_app = Move(cd).(strcat('dir',num2str(ndir)));
            ET_time = NaN(lenM(cd,ndir),2);
            ET_risp = cell(lenM(cd,ndir),1);
            ET_max = NaN(lenM(cd,ndir),1);
            for k = 1:lenM(cd,ndir)
                time_app = M_app(k).(Jtime);
                zero_ind = find(time_app>=0,1,'first');
                t_in = time_app(zero_ind)-abs(tStart);
                tStart_ind =  find(time_app-t_in>=0,1,'first');
                ET_time(k,:) = [time_app(zero_ind)+M_app(k).ET,time_app(end)];
                if ET_time(k,1)<ET_time(k,2)
                    ET_risp{k} = 'SI';
                else
                    ET_risp{k} = 'NO';
                    ET_max(k) = ET_time(k,2);
                end
            end
            ET_session(cd).(strcat('dir',num2str(ndir))) = ET_time;
            ET_rispSess(cd).(strcat('dir',num2str(ndir))) = ET_risp;
            ET_maxSess(cd).(strcat('dir',num2str(ndir))) = ET_max;
        end
    end
    ET_subj(idsess).Sub = ET_session;
    ET_subRis(idsess).Sub = ET_rispSess;
    ET_subMax(idsess).Sub = ET_maxSess;
end

allValues = [];
totalvalues = [];
for i = 1:numel(ET_subMax) 
    subStruct = ET_subMax(i).Sub; 
    for k = 1:length(subStruct)
        subdir = subStruct(k);
        fieldNames = fieldnames(subdir);
        for j = 1:numel(fieldNames)
            values = subdir.(fieldNames{j});
            if isnumeric(values)
                allValues = [allValues; values(~isnan(values))];
            end
        end
        totalvalues = [totalvalues;allValues];
    end
end

% massimo tra tutti i valori non NaN
if ~isempty(totalvalues)
    maxValue = max(totalvalues);
    fprintf('Il valore massimo tra tutti i valori non NaN è: %.4f\n', maxValue);
else
    fprintf('Non ci sono valori validi (non NaN).\n');
end