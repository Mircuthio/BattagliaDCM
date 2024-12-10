%% function TEST_BATTAGLIA_RT_Statistics
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
%% Reaction Time
RTfieldname = 'RT';
RT_S = struct();
RT_Smean = struct();
RT_K = struct();
RT_Kmean = struct();
RT_Sstd = struct();
RT_Kstd = struct();
RT_NaN = struct();
for cd=1:num_cond
    for ndir=1:num_dir
        Move_app = Move(cd).(strcat('dir',num2str(ndir)));
        RT_Sdir = NaN(lenM(cd,ndir),1);
        RT_Kdir = NaN(lenM(cd,ndir),1);
        for k = 1:lenM(cd,ndir)
            RT_Sdir(k,1) = Move_app(k).(strcat(RTfieldname,'_S'));
            RT_Kdir(k,1) = Move_app(k).(strcat(RTfieldname,'_K'));
        end
        RT_S(cd).(strcat('dir',num2str(ndir))) = RT_Sdir;
        RT_K(cd).(strcat('dir',num2str(ndir))) = RT_Kdir;
        RT_NaN(cd).(strcat('dir',num2str(ndir))).len = length(RT_Sdir);
        RT_NaN(cd).(strcat('dir',num2str(ndir))).S = [sum(isnan(RT_Sdir))];
        RT_NaN(cd).(strcat('dir',num2str(ndir))).K = [sum(isnan(RT_Kdir))];
        RT_Smean(cd).(strcat('dir',num2str(ndir))) = mean(RT_Sdir,'omitmissing');
        RT_Kmean(cd).(strcat('dir',num2str(ndir))) = mean(RT_Kdir,'omitmissing');
        RT_Sstd(cd).(strcat('dir',num2str(ndir))) = std(RT_Sdir,'omitmissing');
        RT_Kstd(cd).(strcat('dir',num2str(ndir))) = std(RT_Kdir,'omitmissing');
    end
end

dir_path = strcat('D:\main_scriptDCM\',session_name,'\rng',num2str(par.irng),'\');

%% TEST statistici per DIREZIONE
RT_statistic(1) = RT_S(1);
RT_statistic(2) = RT_K(2);
RT_statistic(3) = RT_S(3);
RT_statistic(4) = RT_K(3);

%% Analisi di Normalità dei dati
skew_QQ = struct();
kurt_QQ = struct();
title_QQ =[condition_name; condition_name{3}];
for cd=1:num_cond+1
    figure;
    for ndir=1:num_dir
        RT_qq = 1000*RT_statistic(cd).(strcat('dir',num2str(ndir)));
        % Test Quantile-Quantile
        subplot(1,num_dir,ndir)
        qqplot(RT_qq);
        title(strcat('Q-Q Plot - dir',num2str(ndir)));
        hold on
        % Asimmetria e Kurtosi
        skew_QQ(cd).(strcat('dir',num2str(ndir))) = skewness(RT_qq);
        kurt_QQ(cd).(strcat('dir',num2str(ndir))) = kurtosis(RT_qq);
        text_string = {strcat('Skewness:',num2str(skewness(RT_qq))),strcat('Kurtosis:',num2str(kurtosis(RT_qq)))};
        x = xlim;
        y = ylim;
        text(x(1),y(2),text_string,'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
    end
    title(title_QQ{cd})
end
%% TEST di Friedman sui 4 gruppi
FriedTest = struct();
for ndir=1:num_dir
    RT_dirFried = {RT_statistic.(strcat('dir',num2str(ndir)))};
    RT_dirFriedclean = cellnan(RT_dirFried);
    RT_dFr = omogdata(RT_dirFriedclean);
    [FriedTest(ndir).p, FriedTest(ndir).tbl, FriedTest(ndir).stats] = friedman(RT_dFr, 1,'off');
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
RT_heat_path = strcat(dir_path,'\RT_Statistic\4Group');

par.savePlotEpsPdfMat.dir_png = strcat(RT_heat_path,'\PNGs\');
par.savePlotEpsPdfMat.dir_pdf = strcat(RT_heat_path,'\PDFs\');
par.savePlotEpsPdfMat.dir_mat = strcat(RT_heat_path,'\MATfiles\');

name_fig = strcat('Heatmaps_4Group_','dir',num2str(ndir));
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
%% TEST di Mann-Whitney & Wilcoxon
alpha = 0.05;
p_value = struct();
for ndir=1:num_dir
    RT_dir = {RT_statistic.(strcat('dir',num2str(ndir)))};
    RT_dirclean = cellnan(RT_dir);
    numGroups = numel(RT_dirclean);
    p2group = nan(numGroups);
    test_name = cell(numGroups);
    for i = 1:numGroups-1
        for j = i+1:numGroups
            if (i == 1 && (j == 2 || j == 4))
                [p2group(i,j), ~] = ranksum(RT_dirclean{i}, RT_dirclean{j},'alpha',alpha);
                test_name{i,j} = 'Mann-Whitney';
            elseif (j == 4 && (i == 2 || i == 3))
                RT_group1 = RT_dirclean{i};
                RT_group2 = RT_dirclean{j};
                if length(RT_group1) ~= length(RT_group2)
                    if length(RT_group1) > length(RT_group2)
                        ndiff = length(RT_group1) - length(RT_group2);
                        idx = randperm(length(RT_group1),ndiff);
                        RT_group1(idx) = [];
                    else
                        ndiff = length(RT_group2) - length(RT_group1);
                        idx = randperm(length(RT_group2),ndiff);
                        RT_group2(idx) = [];
                    end
                end
                [p2group(i,j), ~] = signrank(RT_group1, RT_group2,'alpha',alpha);
                test_name{i,j} = 'Wilcoxon';
            elseif i == 1 && j == 3
                RT_group1 = RT_dirclean{i};
                RT_group2 = RT_dirclean{j};
                if length(RT_group1) ~= length(RT_group2)
                    if length(RT_group1) > length(RT_group2)
                        ndiff = length(RT_group1) - length(RT_group2);
                        idx = randperm(length(RT_group1),ndiff);
                        RT_group1(idx) = [];
                    else
                        ndiff = length(RT_group2) - length(RT_group1);
                        idx = randperm(length(RT_group2),ndiff);
                        RT_group2(idx) = [];
                    end
                end
                [p2group(i,j), ~] = signrank(RT_group1, RT_group2,'alpha',alpha);
                test_name{i,j} = 'Wilcoxon';
            elseif i == 2 && j == 3
                [p2group(i,j), ~] = ranksum(RT_dirclean{i}, RT_dirclean{j},'alpha',alpha);
                test_name{i,j} = 'Mann-Whitney';
            end
        end
    end
    % [pMN12, hMN12] = ranksum(RT_dirclean{1}, RT_dirclean{2});
    % [pMN14, hMN14] = ranksum(RT_dirclean{1}, RT_dirclean{4});
    % [pMN23, hMN23] = ranksum(RT_dirclean{2}, RT_dirclean{3});
    %
    % [pW13, hW13] = signrank(RT_dirclean{1}, RT_dirclean{3});
    % [pW24, hW24] = signrank(RT_dirclean{2}, RT_dirclean{4});
    % [pW34, hW34] = signrank(RT_dirclean{3}, RT_dirclean{4});
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

RT_heat_path = strcat(dir_path,'\RT_Statistic\2Group');

par.savePlotEpsPdfMat.dir_png = strcat(RT_heat_path,'\PNGs\');
par.savePlotEpsPdfMat.dir_pdf = strcat(RT_heat_path,'\PDFs\');
par.savePlotEpsPdfMat.dir_mat = strcat(RT_heat_path,'\MATfiles\');

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
    RT_heat_path = strcat(dir_path,'\RT_Statistic\2Group');

    par.savePlotEpsPdfMat.dir_png = strcat(RT_heat_path,'\PNGs\');
    par.savePlotEpsPdfMat.dir_pdf = strcat(RT_heat_path,'\PDFs\');
    par.savePlotEpsPdfMat.dir_mat = strcat(RT_heat_path,'\MATfiles\');

    name_fig = strcat('Heatmaps_2Group_','dir',num2str(ndir));
    par.savePlotEpsPdfMat.file_name = name_fig;
    savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
end
close all;