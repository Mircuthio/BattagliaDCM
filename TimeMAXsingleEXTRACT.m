
function [TimeMAX_S,TimeMAX_K,TimeMAX_plotMean,TimeMAX_plotStd,max_y,TimeMAX_name] = TimeMAXsingleEXTRACT(data_trials,par)

num_cond = par.num_cond;
num_dir = par.num_dir;
windowSize = par.windowSize;
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
RT = struct();
TimeVelocity = struct();
M_name = struct();
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

TS_smoothSG_Max = struct();
TK_smoothSG_Max = struct();

for cd=1:num_cond
    for ndir = 1:num_dir
        time_app = TimeVelocity(cd).(strcat('dir', num2str(ndir)));
        Vs = MS_smoothSG(cd).(strcat('dir',num2str(ndir)));
        Vk = MK_smoothSG(cd).(strcat('dir',num2str(ndir)));
        Vs_max = NaN(lenM(cd,ndir),1);
        Vk_max = NaN(lenM(cd,ndir),1);
        Ts_max = NaN(lenM(cd,ndir),1);
        Tk_max = NaN(lenM(cd,ndir),1);
        for k=1:lenM(cd,ndir)
            time_appK = time_app(k).new;
            MovOnset_ind = find(time_app(k).new>=0,1,'First');
            time_new = time_appK(MovOnset_ind:end);
            Vs_app = Vs(k).Vel;
            Vk_app = Vk(k).Vel;
            if Vs_app~=0
                [Vs_max(k,1),Ts_maxind]= max(Vs_app(MovOnset_ind:end));
                Ts_max(k,1) = time_new(Ts_maxind);
            else
                Vs_max(k,1) = 0;
                Ts_maxind = 0;
                Ts_max(k,1) = 0;
            end
            if Vk_app~=0
                [Vk_max(k,1),Tk_maxind]  = max(Vk_app(MovOnset_ind:end));
                Tk_max(k,1) = time_new(Tk_maxind);
            else
                Vk_max(k,1) = 0;
                Tk_maxind = 0;
                Tk_max(k,1) = 0;
            end
        end
        TS_smoothSG_Max(cd).(strcat('dir',num2str(ndir))) = Ts_max;
        TK_smoothSG_Max(cd).(strcat('dir',num2str(ndir))) = Tk_max;
    end
end

TimeMAX_S = TS_smoothSG_Max;
TimeMAX_K = TK_smoothSG_Max;
TimeMAX_name = M_name;
%% TEST Statistici per DIREZIONE
TimeMax_bar(1) = TS_smoothSG_Max(1);
TimeMax_bar(2) = TK_smoothSG_Max(2);
TimeMax_bar(3) = TS_smoothSG_Max(3);
TimeMax_bar(4) = TK_smoothSG_Max(3);

%% TimeMAX Max
TimeMAX_Mean = struct();
TimeMAX_Std = struct();
for icond = 1:num_cond+1
    for ndir = 1:ndir
        TimeMAX_calc = TimeMax_bar(icond).(strcat('dir',num2str(ndir)));
        TimeMAX_Mean(icond).(strcat('dir',num2str(ndir))) = mean(TimeMAX_calc);
        TimeMAX_Std(icond).(strcat('dir',num2str(ndir))) = std(TimeMAX_calc);
    end
end
TimeMAX_plotMean = (cell2mat(struct2cell(TimeMAX_Mean')))';
TimeMAX_plotStd = (cell2mat(struct2cell(TimeMAX_Std')))';
offset = max(max(TimeMAX_plotMean+TimeMAX_plotStd))/5;
max_y = max(max(TimeMAX_plotMean+TimeMAX_plotStd))+offset;