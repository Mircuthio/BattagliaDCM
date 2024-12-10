%% function TEST_BATTAGLIA_Jvelocity_plot
clear all; close all;
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

smooth                          = 'SG'; % SG is Savitzky-Golay Filter

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

%% Velocity Grasso
behav_JOINT  = 'H';
signal_type  = 'SUA_BehavProcessed';
data_Grasso  = load([session_name behav_JOINT '_' signal_type]);
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

V_grasso = struct();
for cd=1:num_cond
    for ndir=1:num_dir
        M_app = Move(cd).(strcat('dir',num2str(ndir)));     
        % k = randi(lenM(cd,ndir));
        k =  1;
        time_app = M_app(k).(Jtime);
        zero_ind = find(time_app>=0,1,'first');
        t_in = time_app(zero_ind)-abs(tStart);
        tStart_ind =  find(time_app-t_in>=0,1,'first');
        ET_time = time_app(zero_ind)+M_app(k).ET;
        Jduration_ind = find(time_app-abs(ET_time)>0,1,'First');
        New_Time = time_app(tStart_ind:Jduration_ind);
        RT_timeS = New_Time(find(New_Time>=0,1,'first'))+M_app(k).RT_S;
        RT_timeK = New_Time(find(New_Time>=0,1,'first'))+M_app(k).RT_K;
        RT_indS = find(New_Time-abs(RT_timeS)>0,1,'First');
        RT_indK = find(New_Time-abs(RT_timeK)>0,1,'First');
        % Jduration_ind = 1250;
        MS_appX = double(M_app(k).(Jfieldname)(1,tStart_ind:Jduration_ind));
        MS_appY = double(M_app(k).(Jfieldname)(2,tStart_ind:Jduration_ind));
        MK_appX = double(M_app(k).(Jfieldname)(5,tStart_ind:Jduration_ind));
        MK_appY = double(M_app(k).(Jfieldname)(6,tStart_ind:Jduration_ind));
        MS(cd).(strcat('dir',num2str(ndir))) = [MS_appX;MS_appY];
        MK(cd).(strcat('dir',num2str(ndir))) = [MK_appX;MK_appY];
        % ET(cd).(strcat('dir',num2str(ndir))) = M_app(k).ET;
        RT(cd).(strcat('dir',num2str(ndir))).S = RT_indS;
        RT(cd).(strcat('dir',num2str(ndir))).K = RT_indK;
        TimeVelocity(cd).(strcat('dir',num2str(ndir))) = time_app(tStart_ind:Jduration_ind);
        M_name(cd).(strcat('dir',num2str(ndir))).D = ndir;
        M_name(cd).(strcat('dir',num2str(ndir))).T = M_app(k).trialId;
        K(cd).(strcat('dir',num2str(ndir))) = k;
        K_trialID(cd).(strcat('dir',num2str(ndir))) = M_app(k).trialId;
        % for grasso
        velocityGrasso = data_Grasso.Trials(M_app(k).trialId).v;
        V_grasso(cd).(strcat('dir',num2str(ndir))) = [velocityGrasso(1,tStart_ind:Jduration_ind);velocityGrasso(2,tStart_ind:Jduration_ind)];
    end
end

dir_path = strcat('D:\main_scriptDCM\',session_name,'\rng',num2str(par.irng),'\');

% %% gradient application
% for cd=1:num_cond
%     for ndir = 1:num_dir
%         figure
%         MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
%         MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);
% 
%         MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
%         MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);
% 
%         t = TimeVelocity(cd).(strcat('dir', num2str(ndir)));
%         dt = gradient(t);
%         dx_s = gradient(MS_data_x);
%         dy_s = gradient(MS_data_y);
%         dvx_s = dx_s./dt;
%         dvy_s = dy_s./dt;
% 
%         dv_s = sqrt(dvx_s.^2 + dvy_s.^2);
% 
%         dx_k = gradient(MK_data_x);
%         dy_k = gradient(MK_data_y);
%         dvx_k = dx_k./dt;
%         dvy_k = dy_k./dt;
% 
%         dv_k = sqrt(dvx_k.^2 + dvy_k.^2);
% 
%         plot(1000*t, dv_s, 'b');
%         hold on;
%         plot(1000*t, dv_k, 'g');
%         hold on;
%         xtimeRT_S = 1000*t(RT(cd).(strcat('dir', num2str(ndir))).S);
%         xtimeRT_K = 1000*t(RT(cd).(strcat('dir', num2str(ndir))).K);
%         if ~isnan(xtimeRT_S)
%             xline(xtimeRT_S,'--b','RT-S')
%         end
%         if  ~isnan(xtimeRT_K)
%             xline(xtimeRT_K,'--g','RT-K')
%         end
%         labeltitle = strcat('Dir ', num2str(ndir));
%         % title(labeltitle)
%         legend({'Monkey - S','Monkey - K'},'Location','northwest');
%         textString = sprintf("Chamber: %s\nCondition: %d", data_trials(1).Chamber, cd);
%        dirField = strcat('dir', num2str(ndir)); 
%        textString = sprintf("%s\nD%d-T%d", textString, M_name(cd).(dirField).D, M_name(cd).(dirField).T);
%        text(0.95, 0.9, textString, 'Units', 'normalized', ...
%                  'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
%                  'FontSize', 11, 'Color', 'k');
%         xlabel('time [ms]')
%         ylabel('Velocity Profile')    
%         % sgtitle(append('Condition -',' ',condition_name{cd},newline,labeltitle),'Color',color_name{cd})
%         sgtitle(condition_name{cd},'Color',color_name{cd})
% 
%         set(gcf, 'WindowState', 'maximized');
%         VelocityGradient_path = strcat(dir_path,'\VelocityGradient\');
% 
%         par.savePlotEpsPdfMat.dir_png = strcat(VelocityGradient_path,'\PNGs\');
%         par.savePlotEpsPdfMat.dir_pdf = strcat(VelocityGradient_path,'\PDFs\');
%         par.savePlotEpsPdfMat.dir_mat = strcat(VelocityGradient_path,'\MATfiles\');
% 
%         name_fig = strcat(condition_name{cd},'_VelocityGradient','_Dir ', num2str(ndir));
% 
%         par.savePlotEpsPdfMat.file_name = name_fig;
%         savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
%         % if ~exist(VelocityGradient_path,'dir')
%         %     mkdir(VelocityGradient_path)
%         % end
%         % name_fig = fullfile(VelocityGradient_path,strcat(condition_name{cd},'VelocityGradient.png'));
%         % saveas(gcf,name_fig,'png');
%     end
%     close all
% end
% 
% 
% %% diff and Smoothing application
% for cd=1:num_cond
%     for ndir = 1:num_dir
%         figure
% 
%         MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
%         MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);
% 
%         MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
%         MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);
% 
%         t = TimeVelocity(cd).(strcat('dir', num2str(ndir)));
%         dt = diff(t);
% 
%         MS_dx = diff(MS_data_x);
%         MS_dy = diff(MS_data_y);
% 
%         MS_vx = MS_dx./dt;
%         MS_vy = MS_dy./dt;
%         MS_velocity = sqrt(MS_vx.^2 + MS_vy.^2);
% 
%         windowSize = 50;
%         ofilter = 3; % Filter order Savitzky-Golay
%         %% moving average filter
%         MS_smoothMA = movmean(MS_velocity, windowSize);
% 
%         if ~isnan(MS_velocity)
%             %% Savitzky-Golay filter
%             MS_smoothSG = sgolayfilt(MS_velocity, ofilter, windowSize+1);
%         else
%             MS_smoothSG = [];
%         end
% 
%         MK_dx = diff(MK_data_x);
%         MK_dy = diff(MK_data_y);
% 
%         MK_vx = MK_dx./dt;
%         MK_vy = MK_dy./dt;
%         MK_velocity = sqrt(MK_vx.^2 + MK_vy.^2);
% 
%         %% moving average filter
%         MK_smoothMA = movmean(MK_velocity, windowSize);
% 
%         if ~isnan(MK_velocity)
%             %% Savitzky-Golay filter
%             MK_smoothSG = sgolayfilt(MK_velocity, ofilter, windowSize+1);
%         else
%             MK_smoothSG= [];
%         end
% 
%         if strcmp(smooth,'SG')
%             % figure
%             plot(1000*t(1:end-1), MS_smoothSG, 'b');
%             hold on;
%             plot(1000*t(1:end-1), MK_smoothSG, 'g');
%             hold on;
%         else
%             % figure
%             plot(1000*t(1:end-1), MS_smoothMA, 'b');
%             hold on;
%             plot(1000*t(1:end-1), MK_smoothMA, 'g');
%             hold on;
%         end
%         xtimeRT_S = 1000*t(RT(cd).(strcat('dir', num2str(ndir))).S);
%         xtimeRT_K = 1000*t(RT(cd).(strcat('dir', num2str(ndir))).K);
%         if ~isnan(xtimeRT_S)
%             xline(xtimeRT_S,'--b','RT-S')
%         end
%         if  ~isnan(xtimeRT_K)
%             xline(xtimeRT_K,'--g','RT-K')
%         end
%         labeltitle = strcat('Dir ', num2str(ndir)); 
%         legend({'Monkey - S','Monkey - K'},'Location','northwest');
%         textString = sprintf("Chamber: %s\nCondition: %d", data_trials(1).Chamber, cd);
%         dirField = strcat('dir', num2str(ndir));
%         textString = sprintf("%s\nD%d-T%d", textString, M_name(cd).(dirField).D, M_name(cd).(dirField).T);
%         text(0.95, 0.9, textString, 'Units', 'normalized', ...
%             'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
%             'FontSize', 11, 'Color', 'k');
%         xlabel('time [ms]')
%         ylabel('Velocity Profile')
%         % sgtitle(append('Condition -',' ',condition_name{cd},newline,labeltitle),'Color',color_name{cd})
%         sgtitle(append(condition_name{cd}),'Color',color_name{cd})
% 
%         set(gcf, 'WindowState', 'maximized');
% 
%         VelocitySmooth_path = strcat(dir_path,'\VelocitySmooth\');
% 
%         par.savePlotEpsPdfMat.dir_png = strcat(VelocitySmooth_path,'\PNGs\');
%         par.savePlotEpsPdfMat.dir_pdf = strcat(VelocitySmooth_path,'\PDFs\');
%         par.savePlotEpsPdfMat.dir_mat = strcat(VelocitySmooth_path,'\MATfiles\');
% 
%         name_fig = strcat(condition_name{cd},'_VelocitySmooth','_Dir ', num2str(ndir));
% 
%         par.savePlotEpsPdfMat.file_name = name_fig;
%         savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
%         % text(MS_smoothSG(end), t(end), label_text, ...
%         %     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
%         % text(MK_smoothSG(end), t(end), label_text, ...
%         %     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'g');
%     end
%     % sgtitle(append('Condition -',' ',condition_name{cd}),'Color',color_name{cd})
%     %  Velocitydiff_smooth_path = strcat(dir_path,'\Velocitydiff_smooth');
%     % if ~exist(Velocitydiff_smooth_path,'dir')
%     %     mkdir(Velocitydiff_smooth_path)
%     % end
%     % name_fig = fullfile(Velocitydiff_smooth_path,strcat(condition_name{cd},'Velocitydiff_smooth.png'));
%     % saveas(gcf,name_fig,'png');
%     close all
% end
%% Grasso plot
for cd=1:num_cond
    for ndir = 1:num_dir
        V_S = V_grasso(cd).(strcat('dir', num2str(ndir)))(1,:);
        V_K = V_grasso(cd).(strcat('dir', num2str(ndir)))(2,:);
        figure
        t = TimeVelocity(cd).(strcat('dir', num2str(ndir)));
        plot(1000*t, V_S, 'b');
        hold on;
        plot(1000*t, V_K, 'g');
        hold on;
        xtimeRT_S = 1000*t(RT(cd).(strcat('dir', num2str(ndir))).S);
        xtimeRT_K = 1000*t(RT(cd).(strcat('dir', num2str(ndir))).K);
        if ~isnan(xtimeRT_S)
            xline(xtimeRT_S,'--b','RT-S')
        end
        if  ~isnan(xtimeRT_K)
            xline(xtimeRT_K,'--g','RT-K')
        end
        labeltitle = strcat('Dir ', num2str(ndir));
        % title(labeltitle)
        legend({'Monkey - S','Monkey - K'},'Location','northwest');
        textString = sprintf("Chamber: %s\nCondition: %d", data_trials(1).Chamber, cd);
       dirField = strcat('dir', num2str(ndir)); 
       textString = sprintf("%s\nD%d-T%d", textString, M_name(cd).(dirField).D, M_name(cd).(dirField).T);
       text(0.95, 0.9, textString, 'Units', 'normalized', ...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                 'FontSize', 11, 'Color', 'k');
        xlabel('time [ms]')
        ylabel('Velocity Profile')    
        % sgtitle(append('Condition -',' ',condition_name{cd},newline,labeltitle),'Color',color_name{cd})
        sgtitle(condition_name{cd},'Color',color_name{cd})
    
        set(gcf, 'WindowState', 'maximized');
        VelocityGradient_path = strcat(dir_path,'\GrassoVelocity\');

        par.savePlotEpsPdfMat.dir_png = strcat(VelocityGradient_path,'\PNGs\');
        par.savePlotEpsPdfMat.dir_pdf = strcat(VelocityGradient_path,'\PDFs\');
        par.savePlotEpsPdfMat.dir_mat = strcat(VelocityGradient_path,'\MATfiles\');

        name_fig = strcat(condition_name{cd},'_GrassoVelocity','_Dir ', num2str(ndir));

        par.savePlotEpsPdfMat.file_name = name_fig;
        savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
        close all;
    end
end