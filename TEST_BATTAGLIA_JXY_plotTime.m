%% function TEST_BATTAGLIA_JXY_plotTime
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

% fixed parameters
session_name                    = 'SK009';%'SK004';%{'SK001','SK009'};  % session name
idir                            = 1;        % directions -> 1-8
S                               = filesep;
%% Step 0: arrange trials
fprintf('Step 0: arrange trials\n');
par.BattagliaArrangeTrials              = BattagliaArrangeTrialsParams();
% par.BattagliaArrangeTrials.whichmodel   = 7;        % 7 Model for action session. 5 Model for all sessions
par.BattagliaArrangeTrials.isdemo       = 1;        % getSelectionIndexes(1,1:3);
par.BattagliaArrangeTrials.selS         = 1;
par.BattagliaArrangeTrials.selK         = 1;
par.BattagliaArrangeTrials.session_name = session_name;  % which session
data_trials                             = BattagliaArrangeTrials(par.BattagliaArrangeTrials);

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
M_name = struct();
K = struct();
Time = struct();
K_trialID=struct();
for cd=1:num_cond
    for ndir=1:num_dir
        M_app = Move(cd).(strcat('dir',num2str(ndir)));
        % k = randi(lenM(cd,ndir));
        k=1;
        time_app = M_app(k).(Jtime);
        zero_ind = find(time_app>=0,1,'first');
        ET_time = time_app(zero_ind)+M_app(k).ET;
        Jduration_ind = find(time_app-abs(ET_time)>0,1,'First');
        if isempty(Jduration_ind)
            Jduration_ind = length(time_app);
        end
        % Jduration_ind = 1250;
        New_Time = time_app(zero_ind:Jduration_ind);
        MS_appX = M_app(k).(Jfieldname)(1,zero_ind:Jduration_ind);
        MS_appY = M_app(k).(Jfieldname)(2,zero_ind:Jduration_ind);
        MK_appX = M_app(k).(Jfieldname)(5,zero_ind:Jduration_ind);
        MK_appY = M_app(k).(Jfieldname)(6,zero_ind:Jduration_ind);
        MS(cd).(strcat('dir',num2str(ndir))) = [MS_appX;MS_appY];
        MK(cd).(strcat('dir',num2str(ndir))) = [MK_appX;MK_appY];
        ET(cd).(strcat('dir',num2str(ndir))) = M_app(k).ET;
        M_name(cd).(strcat('dir',num2str(ndir))).D = ndir;
        M_name(cd).(strcat('dir',num2str(ndir))).T = M_app(k).trialId;
        K(cd).(strcat('dir',num2str(ndir))) = k;
        K_trialID(cd).(strcat('dir',num2str(ndir))) = M_app(k).trialId;
        Time(cd).(strcat('dir',num2str(ndir))) = New_Time;
    end
end

dir_path = strcat('D:\main_scriptDCM\',session_name,'\rng',num2str(par.irng),'\');

traj_length_MS = struct();
traj_length_MK = struct();
for cd = 1:num_cond
    for ndir = 1:num_dir
        % Estrazione dei dati x e y per MS (Monkey-S)
        MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
        MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);

        % Estrazione dei dati x e y per MK (Monkey-K)
        MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
        MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);

        % Calcolo della lunghezza della traiettoria (Monkey-S)
        dx_MS = diff(MS_data_x);
        dy_MS = diff(MS_data_y);
        traj_length_MS(cd).(strcat('dir', num2str(ndir))) = sum(sqrt(dx_MS.^2 + dy_MS.^2));

        % Calcolo della lunghezza della traiettoria (Monkey-K)
        dx_MK = diff(MK_data_x);
        dy_MK = diff(MK_data_y);
        traj_length_MK(cd).(strcat('dir', num2str(ndir))) = sum(sqrt(dx_MK.^2 + dy_MK.^2));
    end
end
%% plot cursor x,y position with circle indications

par.plotJxyCircle.condition_name   = condition_name;
par.plotJxyCircle.color_name       = color_name;
par.plotJxyCircle.chamber          = data_trials(1).Chamber;
par.plotJxyCircle.session          = session_name; 
par.plotJxyCircle.dir_path         = dir_path;
par.plotJxyCircle.trjleng          = 0;
par.plotJxyCircle.r_dim            = 1.81;

plotJxyCircle(MS,MK,M_name,par.plotJxyCircle);


%% Grafico distanza dal target singolo

par.plotTargetDistance.condition_name   = condition_name;
par.plotTargetDistance.color_name       = color_name;
par.plotTargetDistance.chamber          = data_trials(1).Chamber;
par.plotTargetDistance.dir_path         = dir_path;
par.plotTargetDistance.trjleng          = 0;
par.plotTargetDistance.circlec          = circlecenter(8);

plotTargetDistance(MS,MK,Time,M_name,par.plotTargetDistance)

%% Grafico distanza dal target all
subplotTargetDistance(MS,MK,Time,M_name,par.plotTargetDistance)

