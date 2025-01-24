%% count number of trials for session
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

dir_path = strcat('D:\main_scriptDCM\','AllSessions','\rng',num2str(par.irng),filesep);


Result = struct();
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

    lenM = NaN(num_cond,num_dir);
    for cd=1:num_cond
        for ndir=1:num_dir
            lenM(cd,ndir) = length(data_trials(Condition{cd,ndir}));
        end
    end

    if sum(sum(lenM)) ~= length(data_trials)
        disp('Error')
    end
    Result(idsess).name = session_name{idsess};
    Result(idsess).chamber = data_trials(1).Chamber;
    Result(idsess).num_trial = lenM;
end
n=1;
good_session = struct();
for i=1:length(Result)
    res_sess = Result(i).num_trial;
    if all(res_sess(:,:) >= 8)
        good_session(n).name_session = Result(i).name;
        good_session(n).chamber = Result(i).chamber;
        good_session(n).lenTrials = Result(i).num_trial;
        n=n+1;
    end
end