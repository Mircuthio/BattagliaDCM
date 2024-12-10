%% function TEST_BATTAGLIA_JXY_plot
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
M_name = struct();
K = struct();
K_trialID=struct();
for cd=1:num_cond
    for ndir=1:num_dir
        M_app = Move(cd).(strcat('dir',num2str(ndir)));     
        % k = randi(lenM(cd,ndir));
        k=1;
        time_app = M_app(k).(Jtime);
        zero_ind = find(time_app>=0,1,'first');
        ET_time = time_app(zero_ind)+M_app(k).ET;
        Jduration_ind = find(time_app-abs(ET_time)>0,1,'First');
        % Jduration_ind = 1250;
        MS_appX = M_app(k).(Jfieldname)(1,zero_ind:Jduration_ind);
        MS_appY = M_app(k).(Jfieldname)(2,zero_ind:Jduration_ind);
        MK_appX = M_app(k).(Jfieldname)(5,zero_ind:Jduration_ind);
        MK_appY = M_app(k).(Jfieldname)(6,zero_ind:Jduration_ind);
        MS(cd).(strcat('dir',num2str(ndir))) = [MS_appX;MS_appY];
        MK(cd).(strcat('dir',num2str(ndir))) = [MK_appX;MK_appY];
        ET(cd).(strcat('dir',num2str(ndir))) = M_app(k).ET;
        M_name(cd).(strcat('dir',num2str(ndir))).D = ndir;
        M_name(cd).(strcat('dir',num2str(ndir))).T = M_app(k).trialId;
        K(cd).(strcat('dir',num2str(ndir))) = k;
        K_trialID(cd).(strcat('dir',num2str(ndir))) = M_app(k).trialId;

    end
end

dir_path = strcat('D:\main_scriptDCM\',session_name,'\rng',num2str(par.irng),'\');
%% plot cursor x,y position
for cd=1:num_cond
    figure
    for ndir = 1:num_dir
        MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
        MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);

        MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
        MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);

        plot(MS_data_x, MS_data_y, 'b');
        xlim([-10 10])
        ylim([-10 10])
        hold on;
        plot(MK_data_x, MK_data_y, 'g');
        xlim([-10 10])
        ylim([-10 10])
        hold on;
        label_text = strcat('Dir ', num2str(ndir)); 
        if ~isnan(MS_data_x)
            text(MS_data_x(end), MS_data_y(end), label_text, ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
        end
        if ~isnan(MK_data_x)
            text(MK_data_x(end), MK_data_y(end), label_text, ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', [0.3, 0.6, 0.4]);
        end
            % Creazione di oggetti fittizi per la legenda
            h1 = plot(nan, nan, 'b');
            h2 = plot(nan, nan, 'g');

            % Aggiungo i cerchi
            r_dim = 1.81;
            addcircle(r_dim);
            % Aggiunta della legenda
            legend([h1 h2], {'Monkey - S', 'Monkey - K'}, 'Location', 'northwest');
             % Positionxyname = append('Condition -',' ',condition_name{cd});
             Positionxyname = append(condition_name{cd});
             title(Positionxyname,'Color',color_name{cd})
    end
   textString = sprintf("Chamber: %s\nCondition: %d", data_trials(1).Chamber, cd);
   for idir = 1:8
       dirField = strcat('dir', num2str(idir)); 
       textString = sprintf("%s\nD%d-T%d", textString, M_name(cd).(dirField).D, M_name(cd).(dirField).T);
   end
       text(0.95, 0.3, textString, 'Units', 'normalized', ...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                 'FontSize', 10, 'Color', 'k');
    set(gcf, 'WindowState', 'maximized');
    XY_path = strcat(dir_path,'\XY');

    par.savePlotEpsPdfMat.dir_png = strcat(XY_path,'\PNGs\');
    par.savePlotEpsPdfMat.dir_pdf = strcat(XY_path,'\PDFs\');
    par.savePlotEpsPdfMat.dir_mat = strcat(XY_path,'\MATfiles\');
    
    name_fig = strcat(condition_name{cd},'_XY');
    par.savePlotEpsPdfMat.file_name = name_fig;
    savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)

    % if ~exist(XY_path,'dir')
    %     mkdir(XY_path)
    % end
    % name_fig = fullfile(XY_path,strcat(condition_name{cd},'_XY.png'));
    % saveas(gcf,name_fig,'png');
    close all;
end
