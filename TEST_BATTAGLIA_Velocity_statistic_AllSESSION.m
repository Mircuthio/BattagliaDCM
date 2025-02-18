%% function TEST_BATTAGLIA_Velocity_statistic_AllSESSION.m
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
    Vel_statistic(1) = MS_smoothSG_Max(1);
    Vel_statistic(2) = MK_smoothSG_Max(2);
    Vel_statistic(3) = MS_smoothSG_Max(3);
    Vel_statistic(4) = MK_smoothSG_Max(3);
    Velocity_session(idsess).cond = Vel_statistic;
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
%% Analisi di Normalità dei dati
skew_QQ = struct();
kurt_QQ = struct();
KS = struct();
title_QQ =[condition_name; condition_name{3}];
for cd=1:num_cond+1
    figure;
    for ndir=1:num_dir
        Velocity_qq = Velocity_total(cd).(strcat('dir',num2str(ndir)));
        [h_ks, p_ks] = kstest((Velocity_qq - mean(Velocity_qq)) / std(Velocity_qq));
        KS(cd).H.(strcat('dir',num2str(ndir))) = h_ks;
        KS(cd).p.(strcat('dir',num2str(ndir))) = p_ks;
        % Test Quantile-Quantile
        subplot(1,num_dir,ndir)
        qqplot(Velocity_qq);
        title(strcat('Q-Q Plot - dir',num2str(ndir)));
        hold on
        % Asimmetria e Kurtosi
        skew_QQ(cd).(strcat('dir',num2str(ndir))) = skewness(Velocity_qq);
        kurt_QQ(cd).(strcat('dir',num2str(ndir))) = kurtosis(Velocity_qq);
        text_string = {strcat('Skewness:',num2str(skewness(Velocity_qq))),strcat('Kurtosis:',num2str(kurtosis(Velocity_qq)))};
        x = xlim;
        y = ylim;
        text(x(1),y(2),text_string,'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
    end
    title(title_QQ{cd})
    set(gcf, 'WindowState', 'maximized');

    Velocity_qq_path = strcat(dir_path,'Velocity_Statistic',filesep,'QQ');

    par.savePlotEpsPdfMat.dir_png = strcat(Velocity_qq_path,'',filesep,'PNGs',filesep,'');
    par.savePlotEpsPdfMat.dir_pdf = strcat(Velocity_qq_path,'',filesep,'PDFs',filesep,'');
    par.savePlotEpsPdfMat.dir_mat = strcat(Velocity_qq_path,'',filesep,'MATfiles',filesep,'');

    name_fig = strcat('QQ','cond',num2str(cd));
    par.savePlotEpsPdfMat.file_name = name_fig;
    savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)

    save(fullfile(Velocity_qq_path,sprintf('KS_irng%d',par.irng)), 'KS');
    close all
end
%% TEST di Friedman sui 4 gruppi
FriedTest = struct();
for ndir=1:num_dir
    Velocity_dirFried = {Velocity_total.(strcat('dir',num2str(ndir)))};
    Velocity_dirFriedclean = cellnan(Velocity_dirFried);
    Velocity_dFr = omogdata(Velocity_dirFriedclean);
    [FriedTest(ndir).p, FriedTest(ndir).tbl, FriedTest(ndir).stats] = friedman(Velocity_dFr, 1,'off');
end
%% plot Friedman test results
FriedGroup1 = {['\bf Act\_S \rm - Obs\_K \newline           ' ...
    'VS \newline \rm Obs\_S - \bf Act\_K \rm \newline           ' ...
    'VS \newline \bf Act\_S \rm - Act\_K \newline           ' ...
    'VS \newline \rm Act\_S - \bf Act\_K']};
FriedGroup2 = {'dir1','dir2', 'dir3', 'dir4', 'dir5', 'dir6', 'dir7', 'dir8'};

figure
h = heatmap(FriedGroup2,FriedGroup1,[FriedTest.p], 'CellLabelFormat', '%.4e', 'ColorLimits', [0 1]);
title(strcat('        All Sessions:','\newline  p value Friedman Test'));
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

Velocity_Fried_path = strcat(dir_path,'Velocity_Statistic',filesep,'4Group');

par.savePlotEpsPdfMat.dir_png = strcat(Velocity_Fried_path,'',filesep,'PNGs',filesep,'');
par.savePlotEpsPdfMat.dir_pdf = strcat(Velocity_Fried_path,'',filesep,'PDFs',filesep,'');
par.savePlotEpsPdfMat.dir_mat = strcat(Velocity_Fried_path,'',filesep,'MATfiles',filesep,'');

name_fig = strcat('Fried_4Group_','dir',num2str(ndir));
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)

save(fullfile(Velocity_Fried_path,sprintf('FriedTest_irng%d',par.irng)), 'FriedTest');
%% TEST di Mann-Whitney & Wilcoxon
alpha = 0.05;
p_value = struct();
for ndir=1:num_dir
    Velocity_dir = {Velocity_total.(strcat('dir',num2str(ndir)))};
    Velocity_dirclean = cellnan(Velocity_dir);
    numGroups = numel(Velocity_dirclean);
    p2group = nan(numGroups);
    test_name = cell(numGroups);
    for i = 1:numGroups-1
        for j = i+1:numGroups
            if (i == 1 && (j == 2 || j == 4))
                [p2group(i,j), ~] = ranksum(Velocity_dirclean{i}, Velocity_dirclean{j},'alpha',alpha);
                test_name{i,j} = 'Mann-Whitney';
            elseif (j == 4 && (i == 2 || i == 3))
                Velocity_group1 = Velocity_dirclean{i};
                Velocity_group2 = Velocity_dirclean{j};
                if length(Velocity_group1) ~= length(Velocity_group2)
                    if length(Velocity_group1) > length(Velocity_group2)
                        ndiff = length(Velocity_group1) - length(Velocity_group2);
                        idx = randperm(length(Velocity_group1),ndiff);
                        Velocity_group1(idx) = [];
                    else
                        ndiff = length(Velocity_group2) - length(Velocity_group1);
                        idx = randperm(length(Velocity_group2),ndiff);
                        Velocity_group2(idx) = [];
                    end
                end
                [p2group(i,j), ~] = signrank(Velocity_group1, Velocity_group2,'alpha',alpha);
                test_name{i,j} = 'Wilcoxon';
            elseif i == 1 && j == 3
                Velocity_group1 = Velocity_dirclean{i};
                Velocity_group2 = Velocity_dirclean{j};
                if length(Velocity_group1) ~= length(Velocity_group2)
                    if length(Velocity_group1) > length(Velocity_group2)
                        ndiff = length(Velocity_group1) - length(Velocity_group2);
                        idx = randperm(length(Velocity_group1),ndiff);
                        Velocity_group1(idx) = [];
                    else
                        ndiff = length(Velocity_group2) - length(Velocity_group1);
                        idx = randperm(length(Velocity_group2),ndiff);
                        Velocity_group2(idx) = [];
                    end
                end
                [p2group(i,j), ~] = signrank(Velocity_group1, Velocity_group2,'alpha',alpha);
                test_name{i,j} = 'Wilcoxon';
            elseif i == 2 && j == 3
                [p2group(i,j), ~] = ranksum(Velocity_dirclean{i}, Velocity_dirclean{j},'alpha',alpha);
                test_name{i,j} = 'Mann-Whitney';
            end
        end
    end
    p_value.(strcat('dir',num2str(ndir))) = p2group;
end
orizz_stripes = 0;
xsquare_stripes =1;
trasv_stripes = 0;
%% plot dei p_value
pvalue_plot = struct();
for ndir=1:num_dir
    p_app = p_value.(strcat('dir',num2str(ndir)));
    p_app(:,1) = [];
    p_app(4,:) = [];
    p_app(2,2) = NaN;
    p_app(1,3) = NaN;
    pvalue_plot.(strcat('dir',num2str(ndir))) = p_app;
end
Groups_nameX = {'      Monkey K \newline Obs\_S - Act\_K \newline      (Solo K)',...
    '       Monkey S \newline Act\_S - Act\_K \newline      (Joint S)',...
    '       Monkey K \newline Act\_S - Act\_K \newline      (Joint K)'};
Groups_nameY = {'      Monkey S \newline Act\_S - Obs\_K \newline      (Solo S)',...
    '       Monkey K \newline Obs\_S - Act\_K \newline      (Solo K)',...
    '       Monkey S \newline Act\_S - Act\_K \newline      (Joint S)'};

figure;
for ndir=1:num_dir
    subplot(2,num_dir/2,ndir)
    curr_pvalue = pvalue_plot.(strcat('dir',num2str(ndir)));
    h = imagesc(curr_pvalue);
    colormap(flipud(gray));
    colorbar;
    clim([0, 1]);
    set(gca, 'Color', [1 1 1]);
    hold on;
    for i = 1:size(curr_pvalue, 1)
        for j = 1:size(curr_pvalue, 2)
            if ~isnan(curr_pvalue(i, j))
                if curr_pvalue(i, j) > 0.5
                    text(j, i, num2str(curr_pvalue(i, j),'%.4f'), 'HorizontalAlignment', 'center', ...
                        'Color', 'white', 'FontWeight', 'bold','FontSize',9);
                elseif curr_pvalue(i, j) <= 0.05
                    if curr_pvalue(i, j) < 0.0001
                        text(j, i, num2str(curr_pvalue(i, j),'%.1e'), 'HorizontalAlignment', 'center', ...
                            'Color', 'black', 'FontWeight', 'bold','FontSize',9);
                    else
                        text(j, i, num2str(curr_pvalue(i, j),'%.4f'), 'HorizontalAlignment', 'center', ...
                            'Color', 'black', 'FontWeight', 'bold','FontSize',9);
                    end
                end
            end
        end
    end
    highlight_cells2 = [1 3; 2 1; 2 2; 3 1; 3 2];
    for k = 1:size(highlight_cells2, 1)
        row = highlight_cells2(k, 1);
        col = highlight_cells2(k, 2);
        rectangle('Position', [col-0.5, row-0.5, 1, 1], ... % Posizione [x, y, larghezza, altezza]
            'EdgeColor', 'k', ... % Colore del bordo
            'LineWidth', 1); % Spessore del bordo
        % Calcolare le coordinate del rettangolo
        x_start = col - 0.5;
        x_end = col + 0.5;
        y_start = row - 0.5;
        y_end = row + 0.5;
        if orizz_stripes ==1
            num_stripes = 30; % Numero di strisce orizzontali
            stripe_height = (y_end - y_start) / num_stripes; % Altezza di ogni striscia
            for i = 1:num_stripes
                stripe_color = [1, 1, 1] * mod(i, 2); % Bianco (1) o Nero (0)
                rectangle('Position', [x_start, y_start + (i-1)*stripe_height, x_end - x_start, stripe_height], ...
                    'FaceColor', stripe_color, 'EdgeColor', 'none');
            end
        elseif xsquare_stripes == 1
            line_width = 1;
            margin = 0.1;
            rectangle('Position', [x_start, y_start, x_end - x_start, y_end - y_start], ...
                'EdgeColor', 'k', 'LineWidth', line_width);
            line([x_start, x_end], [y_start, y_end], 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); % Diagonale \
            line([x_start, x_end], [y_end, y_start], 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); % Diagonale /
        elseif trasv_stripes==1
            rectangle('Position', [col-0.5, row-0.5, 1, 1], ...
                'EdgeColor', 'k', ...
                'LineWidth', 1); % Spessore del bordo
            % Disegna linee trasversali
            num_stripes = 10; % Numero di strisce diagonali
            x_start = col - 0.5;
            x_end = col + 0.5;
            y_start = row - 0.5;
            y_end = row + 0.5;
            for i = 0:num_stripes
                % Calcola la posizione della linea diagonale
                x1 = x_start + i * (x_end - x_start) / num_stripes;
                y1 = y_start;
                x2 = x_start;
                y2 = y_start + i * (y_end - y_start) / num_stripes;
                % Disegna una linea diagonale
                line([x1, x2], [y1, y2], 'Color', 'k', 'LineWidth', 1);
                % % Disegna la linea opposta
                % x1 = x_start + i * (x_end - x_start) / num_stripes;
                % y1 = y_end;
                % x2 = x_end;
                % y2 = y_end - i * (y_end - y_start) / num_stripes;
                % line([x1, x2], [y1, y2], 'Color', 'k', 'LineWidth', 1);
            end
        else            
            num_stripes = 50; % Numero di strisce orizzontali
            stripe_height = (y_end - y_start) / num_stripes; % Altezza di ogni striscia
            stripe_width = (x_end - x_start) / num_stripes; % Larghezza di ogni striscia
            for i = 1:num_stripes
                % Colore alternato per ogni striscia (bianco o nero)
                stripe_color = [1, 1, 1] * mod(i, 2); % Bianco (1) o Nero (0)
                % Aggiungi la striscia
                rectangle('Position', [x_start, y_start + (i-1)*stripe_height, x_end - x_start, stripe_height], ...
                    'FaceColor', stripe_color, 'EdgeColor', 'none');
                rectangle('Position', [x_start + (i-1)*stripe_width, y_start, stripe_width, y_end - y_start], ...
                    'FaceColor', stripe_color, 'EdgeColor', 'none');
            end
        end
    end
    rectangles = [1 1; 3 3; 1 2; 2 3];
    colors = {'m', 'm', 'c', 'c'};
    shrink_factor = 0.1; % Fattore di riduzione per spostare i bordi all'interno
    line_width = 3; % Spessore del bordo

    for k = 1:size(rectangles, 1)
        row = rectangles(k, 1);
        col = rectangles(k, 2);
        % Riduzione delle dimensioni del rettangolo
        rectangle('Position', [col-0.5+shrink_factor/2, row-0.5+shrink_factor/2, 1-shrink_factor, 1-shrink_factor], ...
            'EdgeColor', colors{k}, ...
            'LineWidth', line_width);
    end
    xticks(1:size(curr_pvalue, 2));
    yticks(1:size(curr_pvalue, 1));
    axis square;
    xticklabels(Groups_nameX)
    yticklabels(Groups_nameY)

    title(strcat('Maximum Velocity significance \newline',...
        '                         dir',num2str(ndir)));
    set(gca, 'FontSize', 9);
end
sgtitle('All Sessions:','Interpreter','tex')
set(gcf, 'WindowState', 'maximized');

Velocity_2Group_path = strcat(dir_path,'Velocity_Statistic',filesep,'2Group');

par.savePlotEpsPdfMat.dir_png = strcat(Velocity_2Group_path,'',filesep,'PNGs',filesep,'');
par.savePlotEpsPdfMat.dir_pdf = strcat(Velocity_2Group_path,'',filesep,'PDFs',filesep,'');
par.savePlotEpsPdfMat.dir_mat = strcat(Velocity_2Group_path,'',filesep,'MATfiles',filesep,'');

name_fig = strcat('WilMann_2Group');
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)


save(fullfile(Velocity_2Group_path,sprintf('WillMann_irng%d',par.irng)), 'p_value');
close all;

% plot per direzione
for ndir=1:num_dir
    figure
    curr_pvalue = pvalue_plot.(strcat('dir',num2str(ndir)));
    h = imagesc(curr_pvalue);
    colormap(flipud(gray));
    colorbar;
    clim([0, 1]);
    set(gca, 'Color', [1 1 1]);
    hold on;
    for i = 1:size(curr_pvalue, 1)
        for j = 1:size(curr_pvalue, 2)
            if ~isnan(curr_pvalue(i, j))
                if curr_pvalue(i, j) > 0.5
                    text(j, i, num2str(curr_pvalue(i, j),'%.4f'), 'HorizontalAlignment', 'center', ...
                        'Color', 'white', 'FontWeight', 'bold','FontSize',12);
                elseif curr_pvalue(i, j) <= 0.5
                    if curr_pvalue(i, j) < 0.0001
                        text(j, i, num2str(curr_pvalue(i, j),'%.1e'), 'HorizontalAlignment', 'center', ...
                            'Color', 'black', 'FontWeight', 'bold','FontSize',12);
                    else
                        text(j, i, num2str(curr_pvalue(i, j),'%.4f'), 'HorizontalAlignment', 'center', ...
                            'Color', 'black', 'FontWeight', 'bold','FontSize',12);
                    end
                end
            end
        end
    end
    highlight_cells2 = [1 3; 2 1; 2 2; 3 1; 3 2];
    for k = 1:size(highlight_cells2, 1)
        row = highlight_cells2(k, 1);
        col = highlight_cells2(k, 2);
        rectangle('Position', [col-0.5, row-0.5, 1, 1], ... % Posizione [x, y, larghezza, altezza]
            'EdgeColor', 'k', ... % Colore del bordo
            'LineWidth', 1); % Spessore del bordo
        % Calcolare le coordinate del rettangolo
        x_start = col - 0.5;
        x_end = col + 0.5;
        y_start = row - 0.5;
        y_end = row + 0.5;
        if orizz_stripes ==1
            num_stripes = 30; % Numero di strisce orizzontali
            stripe_height = (y_end - y_start) / num_stripes; % Altezza di ogni striscia
            for i = 1:num_stripes
                stripe_color = [1, 1, 1] * mod(i, 2); % Bianco (1) o Nero (0)
                rectangle('Position', [x_start, y_start + (i-1)*stripe_height, x_end - x_start, stripe_height], ...
                    'FaceColor', stripe_color, 'EdgeColor', 'none');
            end
        elseif xsquare_stripes == 1
            line_width = 1;
            margin = 0.1;
            rectangle('Position', [x_start, y_start, x_end - x_start, y_end - y_start], ...
                'EdgeColor', 'k', 'LineWidth', line_width);
            line([x_start, x_end], [y_start, y_end], 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); % Diagonale \
            line([x_start, x_end], [y_end, y_start], 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); % Diagonale /
        elseif trasv_stripes==1
            rectangle('Position', [col-0.5, row-0.5, 1, 1], ...
                'EdgeColor', 'k', ...
                'LineWidth', 1); % Spessore del bordo
            % Disegna linee trasversali
            num_stripes = 10; % Numero di strisce diagonali
            x_start = col - 0.5;
            x_end = col + 0.5;
            y_start = row - 0.5;
            y_end = row + 0.5;
            for i = 0:num_stripes
                % Calcola la posizione della linea diagonale
                x1 = x_start + i * (x_end - x_start) / num_stripes;
                y1 = y_start;
                x2 = x_start;
                y2 = y_start + i * (y_end - y_start) / num_stripes;
                % Disegna una linea diagonale
                line([x1, x2], [y1, y2], 'Color', 'k', 'LineWidth', 1);
                % % Disegna la linea opposta
                % x1 = x_start + i * (x_end - x_start) / num_stripes;
                % y1 = y_end;
                % x2 = x_end;
                % y2 = y_end - i * (y_end - y_start) / num_stripes;
                % line([x1, x2], [y1, y2], 'Color', 'k', 'LineWidth', 1);
            end
        else            
            num_stripes = 50; % Numero di strisce orizzontali
            stripe_height = (y_end - y_start) / num_stripes; % Altezza di ogni striscia
            stripe_width = (x_end - x_start) / num_stripes; % Larghezza di ogni striscia
            for i = 1:num_stripes
                % Colore alternato per ogni striscia (bianco o nero)
                stripe_color = [1, 1, 1] * mod(i, 2); % Bianco (1) o Nero (0)
                % Aggiungi la striscia
                rectangle('Position', [x_start, y_start + (i-1)*stripe_height, x_end - x_start, stripe_height], ...
                    'FaceColor', stripe_color, 'EdgeColor', 'none');
                rectangle('Position', [x_start + (i-1)*stripe_width, y_start, stripe_width, y_end - y_start], ...
                    'FaceColor', stripe_color, 'EdgeColor', 'none');
            end
        end
    end
    rectangles = [1 1; 3 3; 1 2; 2 3];
    colors = {'m', 'm', 'c', 'c'};
    shrink_factor = 0.1; % Fattore di riduzione per spostare i bordi all'interno
    line_width = 5; % Spessore del bordo

    for k = 1:size(rectangles, 1)
        row = rectangles(k, 1);
        col = rectangles(k, 2);
        % Riduzione delle dimensioni del rettangolo
        rectangle('Position', [col-0.5+shrink_factor/2, row-0.5+shrink_factor/2, 1-shrink_factor, 1-shrink_factor], ...
            'EdgeColor', colors{k}, ...
            'LineWidth', line_width);
    end
    xticks(1:size(curr_pvalue, 2));
    yticks(1:size(curr_pvalue, 1));
    axis square;
    xticklabels(Groups_nameX)
    yticklabels(Groups_nameY)
    title(strcat('                      All Sessions',...
        '\newline Maximum Velocity significance',' dir',num2str(ndir)));
    set(gca, 'FontSize', 14);
    set(gcf, 'WindowState', 'maximized');


    Velocity_2Group_path = strcat(dir_path,'Velocity_Statistic',filesep,'2Group');

    par.savePlotEpsPdfMat.dir_png = strcat(Velocity_2Group_path,'',filesep,'PNGs',filesep,'');
    par.savePlotEpsPdfMat.dir_pdf = strcat(Velocity_2Group_path,'',filesep,'PDFs',filesep,'');
    par.savePlotEpsPdfMat.dir_mat = strcat(Velocity_2Group_path,'',filesep,'MATfiles',filesep,'');

    name_fig = strcat('WillMann_2Group_','dir',num2str(ndir));
    par.savePlotEpsPdfMat.file_name = name_fig;
    savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
end
close all;