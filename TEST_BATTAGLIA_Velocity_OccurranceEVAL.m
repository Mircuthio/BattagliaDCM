%% TEST_BATTAGLIA_Velocity_OccurranceEVAL.m
clear; close all;
par.irng = 10;
rng(par.irng)

windowSize = 50;

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

Velocity_session = struct();
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
    MS = struct();
    MK = struct();
    ET = struct();
    RT = struct();
    TimeVelocity = struct();
    M_name = struct();
    K = struct();
    K_trialID = struct();
    tStart = -0.2;
    for cd=1:num_cond
        for ndir=1:num_dir
            M_app = Move(cd).(strcat('dir',num2str(ndir)));
            % k = randi(lenM(cd,ndir));
            % k =  1;
            RT_ind = struct();
            MS_app = struct();
            MK_app = struct();
            timeSK = struct();
            trialID = NaN(lenM(cd,ndir));
            for k = 1:lenM(cd,ndir)
                time_app = M_app(k).(Jtime);
                zero_ind = find(time_app>=0,1,'first');
                t_in = time_app(zero_ind)-abs(tStart);
                tStart_ind =  find(time_app-t_in>=0,1,'first');
                ET_time = time_app(zero_ind)+M_app(k).ET;
                Jduration_ind = find(time_app-abs(ET_time)>0,1,'First');
                New_Time = time_app(tStart_ind:Jduration_ind);
                RT_timeS = New_Time(find(New_Time>=0,1,'first'))+M_app(k).RT_S;
                RT_timeK = New_Time(find(New_Time>=0,1,'first'))+M_app(k).RT_K;
                RT_ind(k).S = find(New_Time-abs(RT_timeS)>0,1,'First');
                RT_ind(k).K = find(New_Time-abs(RT_timeK)>0,1,'First');
                % Jduration_ind = 1250;
                MS_appX = double(M_app(k).(Jfieldname)(1,tStart_ind:Jduration_ind));
                MS_appY = double(M_app(k).(Jfieldname)(2,tStart_ind:Jduration_ind));
                MK_appX = double(M_app(k).(Jfieldname)(5,tStart_ind:Jduration_ind));
                MK_appY = double(M_app(k).(Jfieldname)(6,tStart_ind:Jduration_ind));
                MS_app(k).XY = [MS_appX;MS_appY];
                MK_app(k).XY = [MK_appX;MK_appY];
                timeSK(k).new = time_app(tStart_ind:Jduration_ind);
                trialID(k,1) = M_app(k).trialId;
                % trialK(k,1) = k;
            end
            MS(cd).(strcat('dir',num2str(ndir))) = MS_app;
            MK(cd).(strcat('dir',num2str(ndir))) = MK_app;
            % ET(cd).(strcat('dir',num2str(ndir))) = M_app(k).ET;
            RT(cd).(strcat('dir',num2str(ndir))).S = RT_ind;
            RT(cd).(strcat('dir',num2str(ndir))).K = RT_ind;
            TimeVelocity(cd).(strcat('dir',num2str(ndir))) = timeSK;
            M_name(cd).(strcat('dir',num2str(ndir))).D = ndir;
            M_name(cd).(strcat('dir',num2str(ndir))).T = trialID;
            % K(cd).(strcat('dir',num2str(ndir))) = k;
            % K_trialID(cd).(strcat('dir',num2str(ndir))) = trialK;
        end
    end

    MS_smoothMA = struct();
    MS_smoothSG = struct();
    MK_smoothMA = struct();
    MK_smoothSG = struct();
    for cd=1:num_cond
        for ndir = 1:num_dir
            MS_data = MS(cd).(strcat('dir', num2str(ndir)));

            MK_data = MK(cd).(strcat('dir', num2str(ndir)));

            timeVel = TimeVelocity(cd).(strcat('dir', num2str(ndir)));
            MS_MA = struct();
            MS_SG = struct();
            MK_MA = struct();
            MK_SG = struct();
            for k=1:lenM(cd,ndir)
                MS_data_x = MS_data(k).XY(1,:);
                MS_data_y = MS_data(k).XY(2,:);
                MK_data_x = MK_data(k).XY(1,:);
                MK_data_y = MK_data(k).XY(2,:);

                dt = diff(timeVel(k).new);

                MS_dx = diff(MS_data_x);
                MS_dy = diff(MS_data_y);

                MS_vx = MS_dx./dt;
                MS_vy = MS_dy./dt;
                MS_velocity = sqrt(MS_vx.^2 + MS_vy.^2);


                ofilter = 3; % Filter order Savitzky-Golay
                %% moving average filter
                MS_MA(k).Vel = movmean(MS_velocity, windowSize);

                if ~isnan(MS_velocity)
                    %% Savitzky-Golay filter
                    MS_SG(k).Vel = sgolayfilt(MS_velocity, ofilter, windowSize+1);
                else
                    MS_SG(k).Vel = 0;
                end

                MK_dx = diff(MK_data_x);
                MK_dy = diff(MK_data_y);

                MK_vx = MK_dx./dt;
                MK_vy = MK_dy./dt;
                MK_velocity = sqrt(MK_vx.^2 + MK_vy.^2);

                %% moving average filter
                MK_MA(k).Vel = movmean(MK_velocity, windowSize);

                if ~isnan(MK_velocity)
                    %% Savitzky-Golay filter
                    MK_SG(k).Vel  = sgolayfilt(MK_velocity, ofilter, windowSize+1);
                else
                    MK_SG(k).Vel = 0;
                end
                MS_smoothMA(cd).(strcat('dir',num2str(ndir))) = MS_MA;
                MS_smoothSG(cd).(strcat('dir',num2str(ndir))) = MS_SG;
                MK_smoothMA(cd).(strcat('dir',num2str(ndir))) = MK_MA;
                MK_smoothSG(cd).(strcat('dir',num2str(ndir))) = MK_SG;
            end
        end
    end

    MS_smoothSG_Max = struct();
    MK_smoothSG_Max = struct();

    for cd=1:num_cond
        for ndir = 1:num_dir
            time_app = TimeVelocity(cd).(strcat('dir', num2str(ndir)));
            Vs = MS_smoothSG(cd).(strcat('dir',num2str(ndir)));
            Vk = MK_smoothSG(cd).(strcat('dir',num2str(ndir)));
            Vs_max = NaN(lenM(cd,ndir),1);
            Vk_max = NaN(lenM(cd,ndir),1);
            for k=1:lenM(cd,ndir)
                MovOnset_ind = find(time_app(k).new>=0,1,'First');
                Vs_app = Vs(k).Vel;
                Vk_app = Vk(k).Vel;
                if Vs_app~=0
                    Vs_max(k,1)= max(Vs_app(MovOnset_ind:end));
                else
                    Vs_max(k,1) = 0;
                end
                if Vk_app~=0
                    Vk_max(k,1)  = max(Vk_app(MovOnset_ind:end));
                else
                    Vk_max(k,1) = 0;
                end
            end
            MS_smoothSG_Max(cd).(strcat('dir',num2str(ndir))) = Vs_max;
            MK_smoothSG_Max(cd).(strcat('dir',num2str(ndir))) = Vk_max;
        end
    end
    %% TEST Statistici per DIREZIONE
    Vel_bar(1) = MS_smoothSG_Max(1);
    Vel_bar(2) = MK_smoothSG_Max(2);
    Vel_bar(3) = MS_smoothSG_Max(3);
    Vel_bar(4) = MK_smoothSG_Max(3);
    Velocity_session(idsess).cond = Vel_bar;
    lenTot(idsess).len = lenM;
end

Velocity_total = struct();
nSession = numel(Velocity_session);
for condIdx = 1:num_cond+1
    for dirIdx = 1:num_dir
        dirName = sprintf('dir%d', dirIdx);
        tempData = struct();
        for nSess = 1:nSession
            if isfield(Velocity_session(nSess).cond(condIdx), dirName)
                tempData(nSess).(dirName) = Velocity_session(nSess).cond(condIdx).(dirName);
            end
        end
        Velocity_total(condIdx).(dirName) = vertcat(tempData.(dirName));
    end
end

par.countVELoccurrance.perc = 50;
[VEL_resultACT,VEL_resultJOINT,VEL_countJOINT,VEL_percJOINT]=countVELoccurrance(Velocity_total,par.countVELoccurrance);

