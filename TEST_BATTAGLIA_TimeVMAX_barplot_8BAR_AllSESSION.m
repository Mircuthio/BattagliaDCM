%% function TEST_BATTAGLIA_TimeVMAX_barplot_8BAR_AllSESSION.m
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

Time_session = struct();
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
    %% TEST Statistici per DIREZIONE
    Time_statistic(1) = TS_smoothSG_Max(1);
    Time_statistic(2) = TK_smoothSG_Max(2);
    Time_statistic(3) = TS_smoothSG_Max(3);
    Time_statistic(4) = TK_smoothSG_Max(3);
    Time_session(idsess).cond = Time_statistic;
    lenTot(idsess).len = lenM;
end
TimeMax_total = struct();
nSession = numel(Time_session);
for condIdx = 1:num_cond+1
    for dirIdx = 1:num_dir
        % Crea il nome del campo per la direzione corrente
        dirName = sprintf('dir%d', dirIdx);
        tempData = struct();
        for nSess = 1:nSession
            if isfield(Time_session(nSess).cond(condIdx), dirName)
                tempData(nSess).(dirName) = Time_session(nSess).cond(condIdx).(dirName);
            end
        end
        % Ora concatenare tutte le matrici nella cella e assegnarle alla struct B
        TimeMax_total(condIdx).(dirName) = vertcat(tempData.(dirName));
    end
end
%% Time Max barplot
categories= {'dir1','dir2', 'dir3', 'dir4', 'dir5', 'dir6', 'dir7', 'dir8'};
TimeMax_Mean = struct();
TimeMax_Std = struct();
for icond = 1:num_cond+1
    for ndir = 1:ndir
        TimeMax_calc = TimeMax_total(icond).(strcat('dir',num2str(ndir)));
        TimeMax_Mean(icond).(strcat('dir',num2str(ndir))) = mean(TimeMax_calc);
        TimeMax_Std(icond).(strcat('dir',num2str(ndir))) = std(TimeMax_calc);
    end
end
TimeMax_plotMean = 1000*(cell2mat(struct2cell(TimeMax_Mean')))';
TimeMax_plotStd = 1000*(cell2mat(struct2cell(TimeMax_Std')))';
offset = max(max(TimeMax_plotMean+TimeMax_plotStd))/5;
max_y = max(max(TimeMax_plotMean+TimeMax_plotStd))+offset;
%% BOXPLOT FIGURE
idcond = 1;

figure
TimeMax_barplot = bar(TimeMax_plotMean');
% infoNaN_plot = infoNaN(cd,:);
% cat_plot = cellfun(@(x, y) [x, '\newline', y], categories, infoNaN_plot, 'UniformOutput', false);
set(gca, 'XTickLabel', categories);
set(gca, 'TickLabelInterpreter', 'tex');
set(gca,'Fontsize',12)
% legend({'Act_S-Obs_K', 'Obs_S-Act_K','Act_S-Act_K','Act_S-Act_K'}, 'Location', 'northwest');
legend({'       Monkey S \newline Act\_S - Obs\_K (Solo S)',...
    '       Monkey K \newline Obs\_S - Act\_K (Solo K)',...
    '       Monkey S \newline Act\_S - Act\_K (Joint S)',...
    '       Monkey K \newline Act\_S - Act\_K (Joint K)'}, ...
    'Location', 'northwest', 'Interpreter', 'tex');
% set(gca, 'XTickLabelRotation', 45);

% orange triple
o = {[1 0.5 0],[1 0.7 0.5],[0.8 0.3 0],[1.0, 0.6, 0.0]};
y = {[1 1 0],[1 1 0.2],[0.8 0.8 0],[1.0, 0.85, 0.0]};
TimeMax_barplot(1).FaceColor = 'b';
TimeMax_barplot(2).FaceColor = 'g';
TimeMax_barplot(3).FaceColor = o{2};
TimeMax_barplot(4).FaceColor = o{1};
ylabel('Time Max [ms]');
ylim([0  max_y+offset]);
TimeMaxtitle_name = append('             All Sessions: ', '\newline', ' Time of maximum velocity'); %,'\newline','*****');
title(TimeMaxtitle_name,'Color','k','Interpreter','tex');
hold on;
x = (1:length(categories));
e1 = errorbar(x - 0.27, TimeMax_plotMean(1,:), TimeMax_plotStd(1,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', 'c', 'LineWidth', 1);
e1.Annotation.LegendInformation.IconDisplayStyle = 'off';

e2 = errorbar(x - 0.09, TimeMax_plotMean(2,:), TimeMax_plotStd(2,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', [0.3, 0.8, 0.6],'LineWidth', 1);
e2.Annotation.LegendInformation.IconDisplayStyle = 'off';

e3 = errorbar(x + 0.09, TimeMax_plotMean(3,:), TimeMax_plotStd(3,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', o{4},'LineWidth', 1);
e3.Annotation.LegendInformation.IconDisplayStyle = 'off';

e4 = errorbar(x + 0.27, TimeMax_plotMean(4,:), TimeMax_plotStd(4,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', o{3},'LineWidth', 1);
e4.Annotation.LegendInformation.IconDisplayStyle = 'off';

% text(0.95, 0.9,strcat('Chamber:',data_trials(1).Chamber), 'Units', 'normalized', ...
%     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
%     'FontSize', 12, 'Color','k');

if idcond==1
    % Testi da aggiungere per ogni gruppo di errorbar
    texts_e1 = repmat({'Solo S'},1,8);
    texts_e2 = repmat({'Solo K'},1,8);
    texts_e3 = repmat({'Joint S'},1,8);
    texts_e4 = repmat({'Joint K'},1,8);
    for i = 1:length(x)
        text(x(i) - 0.27, TimeMax_plotMean(1,i) + TimeMax_plotStd(1,i) + 0.04, texts_e1{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
    for i = 1:length(x)
        text(x(i) - 0.09, TimeMax_plotMean(2,i) + TimeMax_plotStd(2,i) + 0.04, texts_e2{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
    for i = 1:length(x)
        text(x(i) + 0.09, TimeMax_plotMean(3,i) + TimeMax_plotStd(3,i) + 0.04, texts_e3{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
    for i = 1:length(x)
        text(x(i) + 0.27, TimeMax_plotMean(4,i) + TimeMax_plotStd(4,i) + 0.04, texts_e4{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
end
set(gcf, 'WindowState', 'maximized');

TimeMax_barplot_path = strcat(dir_path,'TimeMax_Barplot\8dir');

par.savePlotEpsPdfMat.dir_png = strcat(TimeMax_barplot_path,'\PNGs\');
par.savePlotEpsPdfMat.dir_pdf = strcat(TimeMax_barplot_path,'\PDFs\');
par.savePlotEpsPdfMat.dir_mat = strcat(TimeMax_barplot_path,'\MATfiles\');

name_fig = strcat('TimeMax_8Dir_barplot');
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)

close all;

for ndir=1:num_dir
    condition_name1 = {'      Monkey S',...
        '      Monkey K',...
        '      Monkey S',...
        '      Monkey K'};
    condition_name2 = {'Act\_S - Obs\_K',...
        'Obs\_S - Act\_K',...
        'Act\_S - Act\_K',...
        'Act\_S - Act\_K'};
    condition_name3 = {'      (Solo S)',...
        '      (Solo K)',...
        '      (Joint S)',...
        '      (Joint K)'};


    % orange triple
    o = {[1 0.5 0],[1 0.7 0.5],[0.8 0.3 0],[1.0, 0.6, 0.0]};
    y = {[1 1 0],[1 1 0.2],[0.8 0.8 0],[1.0, 0.85, 0.0]};
    colorplot = {[0 0 1],[0 1 0],o{2},o{1}};
    colorerror = {[0 1 1],[0.3, 0.8, 0.6],o{4},o{3}};

    figure
    TimeMax_pldirMean = TimeMax_plotMean(:,ndir);
    TimeMax_pldirSt = TimeMax_plotStd(:,ndir);
    TimeMax_barplot = bar(TimeMax_pldirMean,'FaceColor', 'flat');
    for k=1:length(TimeMax_pldirMean)
        TimeMax_barplot.CData(k,:) = colorplot{k};
    end
    % infoNaN_plot = infoNaN(cd,:);
    cat_plot = cellfun(@(x, y,z) [x, '\newline', y,'\newline',z], condition_name1, condition_name2,condition_name3, 'UniformOutput', false);
    set(gca, 'XTickLabel', cat_plot);
    set(gca, 'TickLabelInterpreter', 'tex');

    ylabel('Time max [ms]');
    ylim([0  max_y]);
    TimeMaxtitle_name = append('          All Sessions: ',...
        '\newline', 'Time of maximum velocity',...
        '\newline','                 dir ',num2str(ndir));
    title(TimeMaxtitle_name,'Color','k','Interpreter','tex');
    hold on;
    % Testi da aggiungere per ogni gruppo di errorbar
    texts_e = {'Solo S','Solo K','Joint S','Joint K'};
    for k = 1:length(TimeMax_pldirMean)
        x= k;
        e = errorbar(x, TimeMax_pldirMean(k), TimeMax_pldirSt(k), ...
            'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', colorerror{k}, 'LineWidth', 1);
        e.Annotation.LegendInformation.IconDisplayStyle = 'off';
        if idcond==1
            text(x, TimeMax_pldirMean(k) + TimeMax_pldirSt(k) + 5, texts_e{k}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);
        end
    end
    % text(0.95, 0.9,strcat('Chamber:',data_trials(1).Chamber), 'Units', 'normalized', ...
    %     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    %     'FontSize', 12, 'Color','k');
    set(gca, 'FontSize', 14);
    set(gcf, 'WindowState', 'maximized');


    TimeMax_barplot_path = strcat(dir_path,'TimeMax_Barplot\Single_dir');

    par.savePlotEpsPdfMat.dir_png = strcat(TimeMax_barplot_path,'\PNGs\');
    par.savePlotEpsPdfMat.dir_pdf = strcat(TimeMax_barplot_path,'\PDFs\');
    par.savePlotEpsPdfMat.dir_mat = strcat(TimeMax_barplot_path,'\MATfiles\');

    name_fig = strcat('TimeMax_barplot_','dir_',num2str(ndir));
    par.savePlotEpsPdfMat.file_name = name_fig;
    savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
end
close all;


%% TimeMax difference for analysis
TimeMax_Difference = struct();
for ndir=1:num_dir
    TimeMax_diff = TimeMax_plotMean(:,ndir)+TimeMax_plotStd(:,ndir);
    TimeMax_Difference(1).(strcat('dir',num2str(ndir))) = TimeMax_diff(2) - TimeMax_diff(1);
    TimeMax_Difference(2).(strcat('dir',num2str(ndir))) = TimeMax_diff(4) - TimeMax_diff(3);
    TimeMax_Difference(3).(strcat('dir',num2str(ndir))) = TimeMax_diff(3) - TimeMax_diff(1);
    TimeMax_Difference(4).(strcat('dir',num2str(ndir))) = TimeMax_diff(4) - TimeMax_diff(2);
end

TimeMax_DifferenceMEAN = struct();
for ndir=1:num_dir
    TimeMax_diff = TimeMax_plotMean(:,ndir);
    TimeMax_DifferenceMEAN(1).(strcat('dir',num2str(ndir))) = TimeMax_diff(2) - TimeMax_diff(1);
    TimeMax_DifferenceMEAN(2).(strcat('dir',num2str(ndir))) = TimeMax_diff(4) - TimeMax_diff(3);
    TimeMax_DifferenceMEAN(3).(strcat('dir',num2str(ndir))) = TimeMax_diff(3) - TimeMax_diff(1);
    TimeMax_DifferenceMEAN(4).(strcat('dir',num2str(ndir))) = TimeMax_diff(4) - TimeMax_diff(2);
end