%% function TEST_BATTAGLIA_Velocity_statistic.m
clear; close all;
par.irng = 10;
rng(par.irng)


windowSize = 70;

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
%% TEST Statistici per DIREZIONE
Vel_statistic(1) = MS_smoothSG_Max(1);
Vel_statistic(2) = MK_smoothSG_Max(2);
Vel_statistic(3) = MS_smoothSG_Max(3);
Vel_statistic(4) = MK_smoothSG_Max(3);

%% Analisi di Normalità dei dati
skew_QQ = struct();
kurt_QQ = struct();
title_QQ =[condition_name; condition_name{3}];
for cd=1:num_cond+1
    figure;
    for ndir=1:num_dir
        Vel_qq = 1000*Vel_statistic(cd).(strcat('dir',num2str(ndir)));
        % Test Quantile-Quantile
        subplot(1,num_dir,ndir)
        qqplot(Vel_qq);
        title(strcat('Q-Q Plot - dir',num2str(ndir)));
        hold on
        % Asimmetria e Kurtosi
        skew_QQ(cd).(strcat('dir',num2str(ndir))) = skewness(Vel_qq);
        kurt_QQ(cd).(strcat('dir',num2str(ndir))) = kurtosis(Vel_qq);
        text_string = {strcat('Skewness:',num2str(skewness(Vel_qq))),strcat('Kurtosis:',num2str(kurtosis(Vel_qq)))};
        x = xlim;
        y = ylim;
        text(x(1),y(2),text_string,'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
    end
    title(title_QQ{cd})
end

%% TEST di Friedman sui 4 gruppi
FriedTest = struct();
for ndir=1:num_dir
    Vel_dirFried = {Vel_statistic.(strcat('dir',num2str(ndir)))};
    Vel_dirFriedclean = cellnan(Vel_dirFried);
    Vel_dFr = omogdata(Vel_dirFriedclean);
    [FriedTest(ndir).p, FriedTest(ndir).tbl, FriedTest(ndir).stats] = friedman(Vel_dFr, 1,'off');
end
%% plot Friedman test results
FriedGroup1 = {['\bf Act\_S \rm - Obs\_K \newline           ' ...
    'VS \newline \rm Obs\_S - \bf Act\_K \rm \newline           ' ...
    'VS \newline \bf Act\_S \rm - Act\_K \newline           ' ...
    'VS \newline \rm Act\_S - \bf Act\_K']};
FriedGroup2 = {'dir1','dir2', 'dir3', 'dir4', 'dir5', 'dir6', 'dir7', 'dir8'};

figure
h = heatmap(FriedGroup2,FriedGroup1,[FriedTest.p], 'CellLabelFormat', '%.3f', 'ColorLimits', [0 1]);
title(strcat('Session:',session_name,' - Chamber:',data_trials(1).Chamber,'\newline              p value Friedman Test'));
xlabel('Groups');
ylabel('Groups');
h.ColorLimits = [0, 1];
h.MissingDataLabel = '';     % Rimozione etichetta "NaN" per celle bianche
colormap(h, flipud(jet));
h.MissingDataColor = [1 1 1];
% custom_colormap = [linspace(1, 0, 256)', zeros(256, 1), linspace(0, 1, 256)'];
% colormap(h, custom_colormap);
% h.Colormap = custom_colormap;
set(gcf, 'WindowState', 'maximized');
Vel_heat_path = strcat(dir_path,'\Vel_Statistic\4Group\',num2str(windowSize),'wind\');

par.savePlotEpsPdfMat.dir_png = strcat(Vel_heat_path,'\PNGs\');
par.savePlotEpsPdfMat.dir_pdf = strcat(Vel_heat_path,'\PDFs\');
par.savePlotEpsPdfMat.dir_mat = strcat(Vel_heat_path,'\MATfiles\');

name_fig = strcat('Heatmaps_4Group_','dir',num2str(ndir));
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
%% TEST di Mann-Whitney & Wilcoxon
alpha = 0.1;
p_value = struct();
for ndir=1:num_dir
    Vel_dir = {Vel_statistic.(strcat('dir',num2str(ndir)))};
    Vel_dirclean = cellnan(Vel_dir);
    numGroups = numel(Vel_dirclean);
    p2group = nan(numGroups);
    test_name = cell(numGroups);
    for i = 1:numGroups-1
        for j = i+1:numGroups
            if (i == 1 && (j == 2 || j == 4))
                [p2group(i,j), ~] = ranksum(Vel_dirclean{i}, Vel_dirclean{j},'alpha',alpha);
                test_name{i,j} = 'Mann-Whitney';
            elseif (j == 4 && (i == 2 || i == 3))
                Vel_group1 = Vel_dirclean{i};
                Vel_group2 = Vel_dirclean{j};
                if length(Vel_group1) ~= length(Vel_group2)
                    if length(Vel_group1) > length(Vel_group2)
                        ndiff = length(Vel_group1) - length(Vel_group2);
                        idx = randperm(length(Vel_group1),ndiff);
                        Vel_group1(idx) = [];
                    else
                        ndiff = length(Vel_group2) - length(Vel_group1);
                        idx = randperm(length(Vel_group2),ndiff);
                        Vel_group2(idx) = [];
                    end
                end
                [p2group(i,j), ~] = signrank(Vel_group1, Vel_group2,'alpha',alpha);
                test_name{i,j} = 'Wilcoxon';
            elseif i == 1 && j == 3
                Vel_group1 = Vel_dirclean{i};
                Vel_group2 = Vel_dirclean{j};
                if length(Vel_group1) ~= length(Vel_group2)
                    if length(Vel_group1) > length(Vel_group2)
                        ndiff = length(Vel_group1) - length(Vel_group2);
                        idx = randperm(length(Vel_group1),ndiff);
                        Vel_group1(idx) = [];
                    else
                        ndiff = length(Vel_group2) - length(Vel_group1);
                        idx = randperm(length(Vel_group2),ndiff);
                        Vel_group2(idx) = [];
                    end
                end
                [p2group(i,j), ~] = signrank(Vel_group1, Vel_group2,'alpha',alpha);
                test_name{i,j} = 'Wilcoxon';
            elseif i == 2 && j == 3
                [p2group(i,j), ~] = ranksum(Vel_dirclean{i}, Vel_dirclean{j},'alpha',alpha);
                test_name{i,j} = 'Mann-Whitney';
            end
        end
    end
    % [pMN12, hMN12] = ranksum(Vel_dirclean{1}, Vel_dirclean{2});
    % [pMN14, hMN14] = ranksum(Vel_dirclean{1}, Vel_dirclean{4});
    % [pMN23, hMN23] = ranksum(Vel_dirclean{2}, Vel_dirclean{3});
    %
    % [pW13, hW13] = signrank(Vel_dirclean{1}, Vel_dirclean{3});
    % [pW24, hW24] = signrank(Vel_dirclean{2}, Vel_dirclean{4});
    % [pW34, hW34] = signrank(Vel_dirclean{3}, Vel_dirclean{4});
    p_value.(strcat('dir',num2str(ndir))) = p2group;
end


%% plot dei p_value
Groups_name = {'\bf Act\_S \rm - Obs\_K \newline Cond 1', '\rm Obs\_S - \bf Act\_K \newline \rm Cond 2', '\bf Act\_S \rm - Act\_K \newline Cond 3', '\rm Act\_S - \bf Act\_K \newline \rm Cond 3'};
figure;
for ndir=1:num_dir
    subplot(2,num_dir/2,ndir)
    curr_pvalue = p_value.(strcat('dir',num2str(ndir)));
    h = heatmap(Groups_name,Groups_name,curr_pvalue, 'CellLabelFormat', '%.2f', 'ColorLimits', [0 1]);
    title(strcat('p value',' dir',num2str(ndir)));
    xlabel('Groups');
    ylabel('Groups');
    h.ColorLimits = [0, 1]; 
    h.MissingDataLabel = '';     % Rimozione etichetta "NaN" per celle bianche
    colormap(h, flipud(jet));
    h.MissingDataColor = [1 1 1];
    % custom_colormap = [linspace(1, 0, 256)', zeros(256, 1), linspace(0, 1, 256)'];
    % colormap(h, custom_colormap);
    % h.Colormap = custom_colormap;
end
sgtitle(strcat('Session: ',session_name,' - Chamber :',data_trials(1).Chamber),'Interpreter','tex')
set(gcf, 'WindowState', 'maximized');

Vel_heat_path = strcat(dir_path,'\Vel_Statistic\2Group\',num2str(windowSize),'wind\');

par.savePlotEpsPdfMat.dir_png = strcat(Vel_heat_path,'\PNGs\');
par.savePlotEpsPdfMat.dir_pdf = strcat(Vel_heat_path,'\PDFs\');
par.savePlotEpsPdfMat.dir_mat = strcat(Vel_heat_path,'\MATfiles\');

name_fig = strcat('Heatmaps_2Group');
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
close all;
for ndir=1:num_dir
    figure
    curr_pvalue = p_value.(strcat('dir',num2str(ndir)));
    h = heatmap(Groups_name,Groups_name,curr_pvalue, 'CellLabelFormat', '%.2f', 'ColorLimits', [0 1]);
    title(strcat('Session:',session_name,' - Chamber:',data_trials(1).Chamber,'\newline              p value',' dir',num2str(ndir)));
    xlabel('Groups');
    ylabel('Groups');
    h.ColorLimits = [0, 1]; 
    h.MissingDataLabel = '';     % Rimozione etichetta "NaN" per celle bianche
    colormap(h, flipud(jet));
    h.MissingDataColor = [1 1 1];
    % custom_colormap = [linspace(1, 0, 256)', zeros(256, 1), linspace(0, 1, 256)'];
    % colormap(h, custom_colormap);
    % h.Colormap = custom_colormap;
    set(gcf, 'WindowState', 'maximized');
    Vel_heat_path = strcat(dir_path,'\Vel_Statistic\2Group\',num2str(windowSize),'wind\');

    par.savePlotEpsPdfMat.dir_png = strcat(Vel_heat_path,'\PNGs\');
    par.savePlotEpsPdfMat.dir_pdf = strcat(Vel_heat_path,'\PDFs\');
    par.savePlotEpsPdfMat.dir_mat = strcat(Vel_heat_path,'\MATfiles\');

    name_fig = strcat('Heatmaps_2Group_','dir',num2str(ndir));
    par.savePlotEpsPdfMat.file_name = name_fig;
    savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
end
close all;