%% function TEST_BATTAGLIA_TimeVMAX_barplot_8BAR.m
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

% fixed parameters
session_name                    = 'SK004';%'SK004';%{'SK001','SK009'};  % session name
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
    end
end

dir_path = strcat('D:\main_scriptDCM\',session_name,'\rng',num2str(par.irng),'\');

MS_smoothMA = struct();
MS_smoothSG = struct();
MK_smoothMA = struct();
MK_smoothSG = struct();
for cd=1:num_cond
    for ndir = 1:num_dir

        MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
        MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);

        MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
        MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);

        t = TimeVelocity(cd).(strcat('dir', num2str(ndir)));
        dt = diff(t);

        MS_dx = diff(MS_data_x);
        MS_dy = diff(MS_data_y);

        MS_vx = MS_dx./dt;
        MS_vy = MS_dy./dt;
        MS_velocity = sqrt(MS_vx.^2 + MS_vy.^2);


        ofilter = 3; % Filter order Savitzky-Golay
        %% moving average filter
        MS_smoothMA(cd).(strcat('dir',num2str(ndir))) = movmean(MS_velocity, windowSize);

        if ~isnan(MS_velocity)
            %% Savitzky-Golay filter
            MS_smoothSG(cd).(strcat('dir',num2str(ndir))) = sgolayfilt(MS_velocity, ofilter, windowSize+1);
        else
            MS_smoothSG(cd).(strcat('dir',num2str(ndir))) = 0;
        end

        MK_dx = diff(MK_data_x);
        MK_dy = diff(MK_data_y);

        MK_vx = MK_dx./dt;
        MK_vy = MK_dy./dt;
        MK_velocity = sqrt(MK_vx.^2 + MK_vy.^2);

        %% moving average filter
        MK_smoothMA(cd).(strcat('dir',num2str(ndir))) = movmean(MK_velocity, windowSize);

        if ~isnan(MK_velocity)
            %% Savitzky-Golay filter
            MK_smoothSG(cd).(strcat('dir',num2str(ndir))) = sgolayfilt(MK_velocity, ofilter, windowSize+1);
        else
            MK_smoothSG(cd).(strcat('dir',num2str(ndir))) = 0;
        end
    end
end

TimeSG_S = struct();
TimeSG_K = struct();

for cd=1:num_cond
    for ndir = 1:num_dir
        time_app = TimeVelocity(cd).(strcat('dir', num2str(ndir)));
        MovOnset_ind = find(time_app>=0,1,'First');
        Vs = MS_smoothSG(cd).(strcat('dir',num2str(ndir)));
        Vk = MK_smoothSG(cd).(strcat('dir',num2str(ndir)));
        t_new = time_app(MovOnset_ind:end);
        if Vs~=0
            [vs_max,ts] = max(Vs(MovOnset_ind:end));
            TimeSG_S(cd).(strcat('dir',num2str(ndir))) = t_new(ts);
        else
            TimeSG_S(cd).(strcat('dir',num2str(ndir))) = 0;
        end
        if Vk~=0
            [vk_max,tk] = max(Vk(MovOnset_ind:end));
            TimeSG_K(cd).(strcat('dir',num2str(ndir))) = t_new(tk);
        else
            TimeSG_K(cd).(strcat('dir',num2str(ndir))) = 0;
        end
    end
end
%% Time Velocity Max barplot
categories= {'dir1','dir2', 'dir3', 'dir4', 'dir5', 'dir6', 'dir7', 'dir8'};

timeS_SGplotMax = (cell2mat(struct2cell(TimeSG_S')))';
timeK_SGplotMax = (cell2mat(struct2cell(TimeSG_K')))';

%% Vmax plot per 8 direzioni per 4 condizioni
Tmax_plotSG = [timeS_SGplotMax(1,:);timeK_SGplotMax(2,:);timeS_SGplotMax(3,:);timeK_SGplotMax(3,:)];
max_y = max(max(Tmax_plotSG));
%% BOXPLOT FIGURE

idcond = 1;

figure
TmaxSG_barplot = bar(Tmax_plotSG');
% cat_plot = cellfun(@(x, y) [x, '\newline', y], categories, infoNaN_plot, 'UniformOutput', false);
set(gca, 'XTickLabel', categories);
set(gca, 'TickLabelInterpreter', 'tex');
% legend({'Act_S-Obs_K', 'Obs_S-Act_K','Act_S-Act_K','Act_S-Act_K'}, 'Location', 'northwest');
legend({'\bf Act\_S \rm - Obs\_K', '\rm Obs\_S - \bf Act\_K', '\bf Act\_S \rm - Act\_K', '\rm Act\_S - \bf Act\_K'}, ...
    'Location', 'northwest', 'Interpreter', 'tex');

% set(gca, 'XTickLabelRotation', 45);
% orange triple
o = {[1 0.5 0],[1 0.7 0.5],[0.8 0.3 0],[1.0, 0.6, 0.0]};
y = {[1 1 0],[1 1 0.2],[0.8 0.8 0],[1.0, 0.85, 0.0]};
TmaxSG_barplot(1).FaceColor = 'b';
TmaxSG_barplot(2).FaceColor = 'g';
TmaxSG_barplot(3).FaceColor = o{2};
TmaxSG_barplot(4).FaceColor = o{1};
ylabel('T(V_{max}) [s]','Interpreter','tex');
ylim([0  ceil(max_y)+max_y/3]);
TimeMaxtitle_name = append('Time of maximum velocity','\newline','*****');
title(TimeMaxtitle_name,'Color','k');
hold on;
x = (1:length(categories));
% e1 = errorbar(x - 0.27, Vmax_plotSG(1,:), RT_plotStd(1,:), ...
%     'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', 'c', 'LineWidth', 1);
% e1.Annotation.LegendInformation.IconDisplayStyle = 'off';
%
% e2 = errorbar(x - 0.09, Vmax_plotSG(2,:), RT_plotStd(2,:), ...
%     'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', [0.3, 0.8, 0.6],'LineWidth', 1);
% e2.Annotation.LegendInformation.IconDisplayStyle = 'off';
%
% e3 = errorbar(x + 0.09, Vmax_plotSG(3,:), RT_plotStd(3,:), ...
%     'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', o{4},'LineWidth', 1);
% e3.Annotation.LegendInformation.IconDisplayStyle = 'off';
%
% e4 = errorbar(x + 0.27, Vmax_plotSG(4,:), RT_plotStd(4,:), ...
%     'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', o{3},'LineWidth', 1);
% e4.Annotation.LegendInformation.IconDisplayStyle = 'off';

text(0.95, 0.9,strcat('Chamber:',data_trials(1).Chamber), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 12, 'Color','k');
hold off;
set(gcf, 'WindowState', 'maximized');

if idcond==1
        texts_e1 = cell(1,num_dir);
    texts_e2 = cell(1,num_dir);
    texts_e3 = cell(1,num_dir);
    texts_e4 = cell(1,num_dir);
    for ndir=1:num_dir
        dirField = strcat('dir', num2str(ndir));
        texts_e1{1,ndir} = sprintf("%s\nD%d-T%d",'Act S',M_name(1).(dirField).D, M_name(1).(dirField).T);
        texts_e2{1,ndir} = sprintf("%s\nD%d-T%d",'Act K',M_name(2).(dirField).D, M_name(2).(dirField).T);
        texts_e3{1,ndir} = sprintf("%s\nD%d-T%d",'Act S',M_name(3).(dirField).D, M_name(3).(dirField).T);
        texts_e4{1,ndir} = sprintf("%s\nD%d-T%d",'Act K',M_name(3).(dirField).D, M_name(3).(dirField).T);
    end

    % Aggiungi il testo sopra le barre di errore per e1
    for i = 1:length(x)
        text(x(i) - 0.27, Tmax_plotSG(1,i) + 0.05, texts_e1{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end

    % Aggiungi il testo sopra le barre di errore per e2
    for i = 1:length(x)
        text(x(i) - 0.09, Tmax_plotSG(2,i) + 0.05, texts_e2{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end

    % Aggiungi il testo sopra le barre di errore per e3
    for i = 1:length(x)
        text(x(i) + 0.09, Tmax_plotSG(3,i) + 0.05, texts_e3{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end

    % Aggiungi il testo sopra le barre di errore per e4
    for i = 1:length(x)
        text(x(i) + 0.27, Tmax_plotSG(4,i) + 0.05, texts_e4{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
end


Tmax_barplot_path = strcat(dir_path,'\TimeMax\8dir\',num2str(windowSize),'wind\');

par.savePlotEpsPdfMat.dir_png = strcat(Tmax_barplot_path,'\PNGs\');
par.savePlotEpsPdfMat.dir_pdf = strcat(Tmax_barplot_path,'\PDFs\');
par.savePlotEpsPdfMat.dir_mat = strcat(Tmax_barplot_path,'\MATfiles\');

name_fig = strcat(condition_name{cd},'_barplot');
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)


close all