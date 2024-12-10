%% function TEST_BATTAGLIA_RT_barplot_8BAR
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
session_name                    = 'SK011';%'SK004';%{'SK001','SK009'};  % session name 
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
RT_plotStd = 1000*[RT_SplotStd(1,:);RT_KplotStd(2,:);RT_SplotStd(3,:);RT_KplotStd(3,:)];
max_y = max(max(RT_plotMean+RT_plotStd));
%% BOXPLOT FIGURE

idcond = 1;

figure
RT_barplot = bar(RT_plotMean');
% infoNaN_plot = infoNaN(cd,:);
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
RT_barplot(1).FaceColor = 'b';
RT_barplot(2).FaceColor = 'g';
RT_barplot(3).FaceColor = o{2};
RT_barplot(4).FaceColor = o{1};
ylabel('mean RT [ms]');
ylim([0  ceil(max_y/100)*100]);
RTtitle_name = append('Reaction Time','\newline','*****');
title(RTtitle_name,'Color','k');
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


RT_barplot_path = strcat(dir_path,'\RT_Barplot\8dir');

par.savePlotEpsPdfMat.dir_png = strcat(RT_barplot_path,'\PNGs\');
par.savePlotEpsPdfMat.dir_pdf = strcat(RT_barplot_path,'\PDFs\');
par.savePlotEpsPdfMat.dir_mat = strcat(RT_barplot_path,'\MATfiles\');

name_fig = strcat(condition_name{cd},'_barplot');
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)


close all