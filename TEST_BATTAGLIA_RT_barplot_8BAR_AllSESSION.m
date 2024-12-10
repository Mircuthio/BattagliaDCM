%% function TEST_BATTAGLIA_RT_barplot_8BAR_AllSESSION
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


RT_session = struct();
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
    data_trials           = BattagliaArrangeTrials(par.BattagliaArrangeTrials);
    
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
    %% Reaction Time
    RTfieldname = 'RT';
    RT_S = struct();
    RT_K = struct();
    % RT_Smean = struct();
    % RT_Kmean = struct();
    % RT_Sstd = struct();
    % RT_Kstd = struct();
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
        %     RT_Smean(cd).(strcat('dir',num2str(ndir))) = mean(RT_Sdir,'omitmissing');
        %     RT_Kmean(cd).(strcat('dir',num2str(ndir))) = mean(RT_Kdir,'omitmissing');
        %     RT_Sstd(cd).(strcat('dir',num2str(ndir))) = std(RT_Sdir,'omitmissing');
        %     RT_Kstd(cd).(strcat('dir',num2str(ndir))) = std(RT_Kdir,'omitmissing');
        end
    end

    %% Braplot per DIREZIONE
    RT_bar(1) = RT_S(1);
    RT_bar(2) = RT_K(2);
    RT_bar(3) = RT_S(3);
    RT_bar(4) = RT_K(3);
    RT_session(idsess).cond = RT_bar;
    lenTot(idsess).len = lenM;
end
RT_total = struct();
nSession = numel(RT_session); 
for condIdx = 1:num_cond+1
    for dirIdx = 1:num_dir
        dirName = sprintf('dir%d', dirIdx);
        tempData = struct();
        for nSess = 1:nSession
            if isfield(RT_session(nSess).cond(condIdx), dirName)
                tempData(nSess).(dirName) = RT_session(nSess).cond(condIdx).(dirName);
            end
        end
        RT_total(condIdx).(dirName) = vertcat(tempData.(dirName));
    end
end
%% RT barplot
infoNaN = stringNaN(RT_NaN,num_cond,num_dir);
categories= {'dir1','dir2', 'dir3', 'dir4', 'dir5', 'dir6', 'dir7', 'dir8'};
RT_Mean = struct();
RT_Std = struct();
for icond = 1:num_cond+1
    for ndir = 1:ndir
        RT_calc = RT_total(icond).(strcat('dir',num2str(ndir)));
        RT_Mean(icond).(strcat('dir',num2str(ndir))) = 1000*mean(RT_calc);
        RT_Std(icond).(strcat('dir',num2str(ndir))) = 1000*std(RT_calc);
    end
end
RT_plotMean = (cell2mat(struct2cell(RT_Mean')))'; 
RT_plotStd = (cell2mat(struct2cell(RT_Std')))';
max_y = max(max(RT_plotMean+2*RT_plotStd));
%% BOXPLOT FIGURE
idcond = 1;

figure
RT_barplot = bar(RT_plotMean');
% infoNaN_plot = infoNaN(cd,:);
% cat_plot = cellfun(@(x, y) [x, '\newline', y], categories, infoNaN_plot, 'UniformOutput', false);
set(gca, 'XTickLabel', categories);
set(gca, 'TickLabelInterpreter', 'tex');
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
RT_barplot(1).FaceColor = 'b';
RT_barplot(2).FaceColor = 'g';
RT_barplot(3).FaceColor = o{2};
RT_barplot(4).FaceColor = o{1};
ylabel('mean RT [ms]');
ylim([0  ceil(max_y/100)*100]);
RTtitle_name = append('  All Sessions: ', '\newline', ' Reaction Time'); %,'\newline','*****');
title(RTtitle_name,'Color','k','Interpreter','tex');
hold on;
x = (1:length(categories));
e1 = errorbar(x - 0.27, RT_plotMean(1,:), RT_plotStd(1,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', 'c', 'LineWidth', 1);
e1.Annotation.LegendInformation.IconDisplayStyle = 'off';

e2 = errorbar(x - 0.09, RT_plotMean(2,:), RT_plotStd(2,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', [0.3, 0.8, 0.6],'LineWidth', 1);
e2.Annotation.LegendInformation.IconDisplayStyle = 'off';

e3 = errorbar(x + 0.09, RT_plotMean(3,:), RT_plotStd(3,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', o{4},'LineWidth', 1);
e3.Annotation.LegendInformation.IconDisplayStyle = 'off';

e4 = errorbar(x + 0.27, RT_plotMean(4,:), RT_plotStd(4,:), ...
    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', o{3},'LineWidth', 1);
e4.Annotation.LegendInformation.IconDisplayStyle = 'off';

% text(0.95, 0.9,strcat('Chamber:',data_trials(1).Chamber), 'Units', 'normalized', ...
%     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
%     'FontSize', 12, 'Color','k');
hold off;
set(gcf, 'WindowState', 'maximized');

if idcond==1
    % Testi da aggiungere per ogni gruppo di errorbar
    texts_e1 = repmat({'Solo S'},1,8);
    texts_e2 = repmat({'Solo K'},1,8);
    texts_e3 = repmat({'Joint S'},1,8);
    texts_e4 = repmat({'Joint K'},1,8);

    % Aggiungi il testo sopra le barre di errore per e1
    for i = 1:length(x)
        text(x(i) - 0.27, RT_plotMean(1,i) + RT_plotStd(1,i) + 0.05, texts_e1{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end

    % Aggiungi il testo sopra le barre di errore per e2
    for i = 1:length(x)
        text(x(i) - 0.09, RT_plotMean(2,i) + RT_plotStd(2,i) + 0.05, texts_e2{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end

    % Aggiungi il testo sopra le barre di errore per e3
    for i = 1:length(x)
        text(x(i) + 0.09, RT_plotMean(3,i) + RT_plotStd(3,i) + 0.05, texts_e3{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end

    % Aggiungi il testo sopra le barre di errore per e4
    for i = 1:length(x)
        text(x(i) + 0.27, RT_plotMean(4,i) + RT_plotStd(4,i) + 0.05, texts_e4{i}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
end


RT_barplot_path = strcat(dir_path,'RT_Barplot\8dir');

par.savePlotEpsPdfMat.dir_png = strcat(RT_barplot_path,'\PNGs\');
par.savePlotEpsPdfMat.dir_pdf = strcat(RT_barplot_path,'\PDFs\');
par.savePlotEpsPdfMat.dir_mat = strcat(RT_barplot_path,'\MATfiles\');

name_fig = strcat('RT_8Dir_barplot');
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)

close all;

for ndir=1:num_dir
    condition_name1 = {'       Monkey S',...
        '       Monkey K',...
        '       Monkey S',...
        '       Monkey K'};
    condition_name2 = {'Act\_S - Obs\_K',...
        'Obs\_S - Act\_K',...
        'Act\_S - Act\_K',...
        'Act\_S - Act\_K'};
    condition_name3 = {'       (Solo S)',...
        '       (Solo K)',...
        '       (Joint S)',...
        '       (Joint K)'};


    % orange triple
    o = {[1 0.5 0],[1 0.7 0.5],[0.8 0.3 0],[1.0, 0.6, 0.0]};
    y = {[1 1 0],[1 1 0.2],[0.8 0.8 0],[1.0, 0.85, 0.0]};
    colorplot = {[0 0 1],[0 1 0],o{2},o{1}};
    colorerror = {[0 1 1],[0.3, 0.8, 0.6],o{4},o{3}};

    figure
    RT_pldirMean = RT_plotMean(:,ndir);
    RT_pldirSt = RT_plotStd(:,ndir);
    RT_barplot = bar(RT_pldirMean,'FaceColor', 'flat');
    for k=1:length(RT_pldirMean)
        RT_barplot.CData(k,:) = colorplot{k};
    end
    % infoNaN_plot = infoNaN(cd,:);
    cat_plot = cellfun(@(x, y,z) [x, '\newline', y,'\newline',z], condition_name1, condition_name2,condition_name3, 'UniformOutput', false);
    set(gca, 'XTickLabel', cat_plot);
    set(gca, 'TickLabelInterpreter', 'tex');

    ylabel('mean RT [ms]');
    ylim([0  ceil(max_y/100)*100]);
    RTtitle_name = append('  All Sessions: ', '\newline', ' Reaction Time',...
        '\newline','         dir ',num2str(ndir));
    title(RTtitle_name,'Color','k','Interpreter','tex');
    hold on;
    % Testi da aggiungere per ogni gruppo di errorbar
    texts_e = {'Solo S','Solo K','Joint S','Joint K'};
    for k = 1:length(RT_pldirMean)
        x= k;
        e = errorbar(x, RT_pldirMean(k), RT_pldirSt(k), ...
            'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', colorerror{k}, 'LineWidth', 1);
        e.Annotation.LegendInformation.IconDisplayStyle = 'off';
        if idcond==1
            text(x, RT_pldirMean(k) + RT_pldirSt(k) + 5, texts_e{k}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);
        end
    end
    % text(0.95, 0.9,strcat('Chamber:',data_trials(1).Chamber), 'Units', 'normalized', ...
    %     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    %     'FontSize', 12, 'Color','k');
    set(gca, 'FontSize', 14);
    set(gcf, 'WindowState', 'maximized');


    RT_barplot_path = strcat(dir_path,'RT_Barplot\Single_dir');

    par.savePlotEpsPdfMat.dir_png = strcat(RT_barplot_path,'\PNGs\');
    par.savePlotEpsPdfMat.dir_pdf = strcat(RT_barplot_path,'\PDFs\');
    par.savePlotEpsPdfMat.dir_mat = strcat(RT_barplot_path,'\MATfiles\');

    name_fig = strcat('RT_barplot_','dir_',num2str(ndir));
    par.savePlotEpsPdfMat.file_name = name_fig;
    savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
end
close all;


%% RT difference for analysis
RT_Difference = struct();
for ndir=1:num_dir
    RT_diff = RT_plotMean(:,ndir)+RT_plotStd(:,ndir);
    RT_Difference(1).(strcat('dir',num2str(ndir))) = RT_diff(2) - RT_diff(1);
    RT_Difference(2).(strcat('dir',num2str(ndir))) = RT_diff(4) - RT_diff(3);
    RT_Difference(3).(strcat('dir',num2str(ndir))) = RT_diff(3) - RT_diff(1);
    RT_Difference(4).(strcat('dir',num2str(ndir))) = RT_diff(4) - RT_diff(2);
end