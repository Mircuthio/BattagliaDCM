%% function TEST_BATTAGLIA_Velocity_barplot_MeanBAR.m
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
        trialid = struct();
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
            trialK(k,1) = k;
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

dir_path = strcat('D:\main_scriptDCM\',session_name,'\rng',num2str(par.irng),'\');

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
%% Velocity Max barplot
categories= {'dir1','dir2', 'dir3', 'dir4', 'dir5', 'dir6', 'dir7', 'dir8'};
MS_SGplotMax = struct();
MK_SGplotMax = struct();
MS_SGplotStd = struct();
MK_SGplotStd = struct();
for cd=1:num_cond
    for ndir = 1:num_dir
        MS_SGplotMax(cd).(strcat('dir',num2str(ndir))) = mean(MS_smoothSG_Max(cd).(strcat('dir',num2str(ndir))));
        MK_SGplotMax(cd).(strcat('dir',num2str(ndir))) = mean(MK_smoothSG_Max(cd).(strcat('dir',num2str(ndir))));
        MS_SGplotStd(cd).(strcat('dir',num2str(ndir))) = std(MS_smoothSG_Max(cd).(strcat('dir',num2str(ndir))));
        MK_SGplotStd(cd).(strcat('dir',num2str(ndir))) = std(MK_smoothSG_Max(cd).(strcat('dir',num2str(ndir))));        
    end
end

MS_SGplotMax = (cell2mat(struct2cell(MS_SGplotMax')))';
MK_SGplotMax = (cell2mat(struct2cell(MK_SGplotMax')))';
MS_SGplotStd = (cell2mat(struct2cell(MS_SGplotStd')))';
MK_SGplotStd = (cell2mat(struct2cell(MK_SGplotStd')))';
%% Vmax plot per 8 direzioni per 4 condizioni
Vmax_plotSG = [MS_SGplotMax(1,:);MK_SGplotMax(2,:);MS_SGplotMax(3,:);MK_SGplotMax(3,:)];
Vmax_plotSGstd = [MS_SGplotStd(1,:);MK_SGplotStd(2,:);MS_SGplotStd(3,:);MK_SGplotStd(3,:)];

max_y = max(max(Vmax_plotSG))+max(max(Vmax_plotSGstd));

%% BOXPLOT FIGURE

idcond = 1;

figure
VmaxSG_barplot = bar(Vmax_plotSG');
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
VmaxSG_barplot(1).FaceColor = 'b';
VmaxSG_barplot(2).FaceColor = 'g';
VmaxSG_barplot(3).FaceColor = o{2};
VmaxSG_barplot(4).FaceColor = o{1};
ylabel('Mean of Vpeak [#]');
ylim([0  (ceil(max_y/100)+max_y/300)*100]);
Velocitytitle_name = append('Mean of maximum Velocity','\newline','*****');
title(Velocitytitle_name,'Color','k');
hold on;
x = (1:length(categories));
e1 = errorbar(x - 0.27, Vmax_plotSG(1,:), Vmax_plotSGstd(1,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', 'c', 'LineWidth', 1);
e1.Annotation.LegendInformation.IconDisplayStyle = 'off';

e2 = errorbar(x - 0.09, Vmax_plotSG(2,:), Vmax_plotSGstd(2,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', [0.3, 0.8, 0.6],'LineWidth', 1);
e2.Annotation.LegendInformation.IconDisplayStyle = 'off';

e3 = errorbar(x + 0.09, Vmax_plotSG(3,:), Vmax_plotSGstd(3,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', o{4},'LineWidth', 1);
e3.Annotation.LegendInformation.IconDisplayStyle = 'off';

e4 = errorbar(x + 0.27, Vmax_plotSG(4,:), Vmax_plotSGstd(4,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', o{3},'LineWidth', 1);
e4.Annotation.LegendInformation.IconDisplayStyle = 'off';

text(0.95, 0.9,strcat('Chamber:',data_trials(1).Chamber), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 12, 'Color','k');
hold off;
set(gcf, 'WindowState', 'maximized');

if idcond==1
    % Testi da aggiungere per ogni gruppo di errorbar
    texts_e1 = repmat({'Act S'},1,8);
    texts_e2 = repmat({'Act K'},1,8);
    texts_e3 = repmat({'Act S'},1,8);
    texts_e4 = repmat({'Act K'},1,8);

    % Aggiungi il testo sopra le barre di errore per e1
    for i = 1:length(x)
        text(x(i) - 0.27, Vmax_plotSG(1,i) + 0.05, texts_e1{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end

    % Aggiungi il testo sopra le barre di errore per e2
    for i = 1:length(x)
        text(x(i) - 0.09, Vmax_plotSG(2,i) + 0.05, texts_e2{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end

    % Aggiungi il testo sopra le barre di errore per e3
    for i = 1:length(x)
        text(x(i) + 0.09, Vmax_plotSG(3,i) + 0.05, texts_e3{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end

    % Aggiungi il testo sopra le barre di errore per e4
    for i = 1:length(x)
        text(x(i) + 0.27, Vmax_plotSG(4,i) + 0.05, texts_e4{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
end


Vmax_barplot_path = strcat(dir_path,'\VelocitySmooth\MeanDir\',num2str(windowSize),'wind\');

par.savePlotEpsPdfMat.dir_png = strcat(Vmax_barplot_path,'\PNGs\');
par.savePlotEpsPdfMat.dir_pdf = strcat(Vmax_barplot_path,'\PDFs\');
par.savePlotEpsPdfMat.dir_mat = strcat(Vmax_barplot_path,'\MATfiles\');

name_fig = strcat(condition_name{cd},'_barplot');
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)


close all