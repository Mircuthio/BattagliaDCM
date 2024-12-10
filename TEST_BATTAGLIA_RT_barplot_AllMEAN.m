%% function TEST_BATTAGLIA_RT_barplot_AllMEAN
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
%% RT barplot
infoNaN = stringNaN(RT_NaN,num_cond,num_dir);
categories= {'dir1','dir2', 'dir3', 'dir4', 'dir5', 'dir6', 'dir7', 'dir8'};

RT_SplotMean = (cell2mat(struct2cell(RT_Smean')))';
RT_KplotMean = (cell2mat(struct2cell(RT_Kmean')))';
RT_SplotStd = (cell2mat(struct2cell(RT_Sstd')))';
RT_KplotStd = (cell2mat(struct2cell(RT_Kstd')))';


%% RT plot per 8 direzioni per 4 condizioni
RT_plotMean = 1000*[RT_SplotMean(1,:);RT_KplotMean(2,:);RT_SplotMean(3,:);RT_KplotMean(3,:)];
% RT_plotStd = 1000*[RT_SplotStd(1,:);RT_KplotStd(2,:);RT_SplotStd(3,:);RT_KplotStd(3,:)];

RT_mean     = mean(RT_plotMean,2);
RT_std      = std(RT_plotMean,0,2);
max_y = max(max(RT_mean+RT_std));
%% BOXPLOT FIGURE
idcond = 1;

condition_name1 = {'Act S';'Obs S';'Act S';'Act K'};
condition_name2 = {'Obs K';'Act K';'Act K';'Act S'};

% orange triple
o = {[1 0.5 0],[1 0.7 0.5],[0.8 0.3 0],[1.0, 0.6, 0.0]};
y = {[1 1 0],[1 1 0.2],[0.8 0.8 0],[1.0, 0.85, 0.0]};
colorplot = {[0 0 1],[0 1 0],o{2},o{1}};
colorerror = {[0 1 1],[0.3, 0.8, 0.6],o{4},o{3}};

figure
RT_barplot = bar(RT_mean,'FaceColor', 'flat');
for k=1:length(RT_mean)
    RT_barplot.CData(k,:) = colorplot{k};
end
% infoNaN_plot = infoNaN(cd,:);
cat_plot = cellfun(@(x, y) [x, '\newline', y], condition_name1, condition_name2, 'UniformOutput', false);
set(gca, 'XTickLabel', cat_plot);
set(gca, 'TickLabelInterpreter', 'tex');
% legend({'Act_S-Obs_K', 'Obs_S-Act_K','Act_S-Act_K','Act_S-Act_K'}, 'Location', 'northwest');
% legend({'\bf Act\_S \rm - Obs\_K', '\rm Obs\_S - \bf Act\_K', '\bf Act\_S \rm - Act\_K', '\rm Act\_S - \bf Act\_K'}, ...
%        'Location', 'northwest', 'Interpreter', 'tex');

% set(gca, 'XTickLabelRotation', 45);

ylabel('mean RT [ms]');
ylim([0  ceil(max_y/100)*100]);
RTtitle_name = append('Mean Reaction Time','\newline','*****');
title(RTtitle_name,'Color','k');
hold on;
% Testi da aggiungere per ogni gruppo di errorbar
texts_e = {'Act S','Act K','Act S','Act K'};
for k = 1:length(RT_mean)
    x= k;
    e = errorbar(x, RT_mean(k), RT_std(k), ...
        'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', colorerror{k}, 'LineWidth', 1);
    e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    if idcond==1
        text(x, RT_mean(k) + RT_std(k) + 5, texts_e{k}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
    end
end
text(0.95, 0.9,strcat('Chamber:',data_trials(1).Chamber), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 12, 'Color','k');
hold off;
set(gcf, 'WindowState', 'maximized');


RT_barplot_path = strcat(dir_path,'\RT_Barplot\AllMean');

par.savePlotEpsPdfMat.dir_png = strcat(RT_barplot_path,'\PNGs\');
par.savePlotEpsPdfMat.dir_pdf = strcat(RT_barplot_path,'\PDFs\');
par.savePlotEpsPdfMat.dir_mat = strcat(RT_barplot_path,'\MATfiles\');

name_fig = strcat(condition_name{cd},'_barplot');
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)


close all