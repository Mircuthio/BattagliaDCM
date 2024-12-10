%% function TEST_BATTAGLIA_DCM_JXYEXY
clear all
rng(10)
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
color_name = {[0 0 0.5],[0 0.5 0],[0.1 0.6 0]};

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
for cd=1:num_cond
    for ndir=1:num_dir
        M_app = Move(cd).(strcat('dir',num2str(ndir)));     
        k = randi(lenM(cd,ndir));
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
dir_path = strcat('D:\main_scriptDCM\',session_name);

%% RT barplot
infoNaN = stringNaN(RT_NaN,num_cond,num_dir);
categories= {'dir1','dir2', 'dir3', 'dir4', 'dir5', 'dir6', 'dir7', 'dir8'};

    RT_SplotMean = (cell2mat(struct2cell(RT_Smean')))';
    RT_KplotMean = (cell2mat(struct2cell(RT_Kmean')))';
    RT_SplotStd = (cell2mat(struct2cell(RT_Sstd')))';
    RT_KplotStd = (cell2mat(struct2cell(RT_Kstd')))';
    maxS = max(max(RT_SplotMean+RT_SplotStd));
    maxK = max(max(RT_KplotMean+RT_KplotStd));

for cd = 1:num_cond
    RT_plotMean = 1000*[RT_SplotMean(cd,:);RT_KplotMean(cd,:)];
    RT_plotStd = 1000*[RT_SplotStd(cd,:);RT_KplotStd(cd,:)];
    figure
    RT_barplot = bar(RT_plotMean');
    infoNaN_plot = infoNaN(cd,:); 
    cat_plot = cellfun(@(x, y) [x, '\newline', y], categories, infoNaN_plot, 'UniformOutput', false);
    set(gca, 'XTickLabel', cat_plot);
    set(gca, 'TickLabelInterpreter', 'tex');
    legend({'Monkey - S', 'Monkey - K'}, 'Location', 'northwest'); 
    % set(gca, 'XTickLabelRotation', 45); 
    RT_barplot(1).FaceColor = 'b'; 
    RT_barplot(2).FaceColor = 'g'; 
    ylabel('RT [ms]');
    ylim([0  round(max(1000*[maxS;maxK])/100)*100]);
    RTtitle_name = append('Reaction Time','\newline',condition_name{cd});
    title(RTtitle_name,'Color',color_name{cd});
    hold on;
    x = (1:length(categories));
    e1 = errorbar(x - 0.15, RT_plotMean(1,:), RT_plotStd(1,:), ...
                  'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', 'c', 'LineWidth', 1);
    e1.Annotation.LegendInformation.IconDisplayStyle = 'off'; 

    e2 = errorbar(x + 0.15, RT_plotMean(2,:), RT_plotStd(2,:), ...
                  'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', [0.3, 0.8, 0.6],'LineWidth', 1);
    e2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold off;
    RT_barplot_path = strcat(dir_path,'\RT_Barplot');
    if ~exist(RT_barplot_path,'dir')
        mkdir(RT_barplot_path)
    end
    name_fig = fullfile(RT_barplot_path,strcat(condition_name{cd},'_barplot.png'));
    saveas(gcf,name_fig,'png');
end

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
        text(MS_data_x(end), MS_data_y(end), label_text, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
        text(MK_data_x(end), MK_data_y(end), label_text, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'g');
            % Creazione di oggetti fittizi per la legenda
            h1 = plot(nan, nan, 'b');
            h2 = plot(nan, nan, 'g');
            r_dim = 1.81;
            addcircle(r_dim);
            % Aggiunta della legenda
            legend([h1 h2], {'Monkey - S', 'Monkey - K'}, 'Location', 'best');
    end
    Positionxyname = append('Condition -',' ',condition_name{cd});
    title(Positionxyname,'Color',color_name{cd})
     XY_path = strcat(dir_path,'\XY');
    if ~exist(XY_path,'dir')
        mkdir(XY_path)
    end
    name_fig = fullfile(XY_path,strcat(condition_name{cd},'_XY.png'));
    saveas(gcf,name_fig,'png');
end

% %% plot cursor x,y position in subplot
% for cd=1:num_cond
%     figure
%     for ndir = 1:num_dir
%         MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
%         MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);
%         MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
%         MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);
%         subplot(2,4,ndir)
%         plot(MS_data_x, MS_data_y, 'b');
%         hold on;
%         plot(MK_data_x, MK_data_y, 'g');
%         hold on;
%         labeltitle = strcat('Dir ', num2str(ndir));
%         title(labeltitle)
%         label_text = strcat('Dir ', num2str(ndir)); 
%         text(MS_data_x(end), MS_data_y(end), label_text, ...
%             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
%         text(MK_data_x(end), MK_data_y(end), label_text, ...
%             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'g');
%             % Creazione di oggetti fittizi per la legenda
%             h1 = plot(nan, nan, 'b');
%             h2 = plot(nan, nan, 'g');
%             r_dim = 1.5;
%             addcircle(r_dim);
%             % Aggiunta della legenda
%             legend([h1 h2], {'Monkey - S', 'Monkey - K',''}, 'Location', 'best');
%     end
%     Positionxyname = append('ConRT_S_datadition -',' ',condition_name{cd});
%     title(Positionxyname,'Color',color_name{cd})
%      XY_subplot_path = strcat(dir_path,'\XY_subplot');
%     if ~exist(XY_subplot_path,'dir')
%         mkdir(XY_subplot_path)
%     end
%     name_fig = fullfile(XY_subplot_path,strcat(condition_name{cd},'XY_subplot.png'));
%     saveas(gcf,name_fig,'png');
% end

%% plot cursor x,y position with RT
if exist('plot_withRT','var')
    for cd=1:num_cond
        figure
        for ndir = 1:num_dir
            % time = ;
            MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
            RT_S_data = RT_Smean(cd).(strcat('dir', num2str(ndir)))(1,:);
            RT_S_xy = find(time-abs(RT_S_data)>0,1,'First');
            MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);

            MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
            RT_K_data = RT_Kmean(cd).(strcat('dir', num2str(ndir)))(1,:);
            RT_K_xy = find(time-abs(RT_K_data)>0,1,'First');
            MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);
            % RT cycle
            plot(MS_data_x(1:RT_S_xy), MS_data_y(1:RT_S_xy), 'r'); % RT mean dir
            hold on;
            plot(MS_data_x(RT_S_xy+1:end), MS_data_y(RT_S_xy+1:end), 'b');
            hold on;
            plot(MK_data_x(1:RT_K_xy), MK_data_y(1:RT_K_xy), 'k'); % RT mean dir
            hold on;
            plot(MK_data_x(RT_K_xy+1:end), MK_data_y(RT_K_xy+1:end), 'g');
            label_text = strcat('Dir ', num2str(ndir));
            text(MS_data_x(end), MS_data_y(end), label_text, ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
            text(MK_data_x(end), MK_data_y(end), label_text, ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'g');
            % Creazione di oggetti fittizi per la legenda
            h1 = plot(nan, nan, 'b');
            h2 = plot(nan, nan, 'g');
            r_dim = 1.5;
            addcircle(r_dim);
            % Aggiunta della legenda
            legend([h1 h2], {'Monkey - S', 'Monkey - K',''}, 'Location', 'best');
        end
        Positionxyname = append('Condition -',' ',condition_name{cd});
        title(Positionxyname,'Color',color_name{cd})
        XY_RT_path = strcat(dir_path,'\XY_RT');
        if ~exist(XY_RT_path,'dir')
            mkdir(XY_RT_path)
        end
        name_fig = fullfile(XY_RT_path,strcat(condition_name{cd},'XY_RT.png'));
        saveas(gcf,name_fig,'png');
    end
end

% %% plot cursor x,y position with RT in subplot 
% time = linspace(0,1,1000);
% for cd=1:num_cond
%     figure
%     for ndir = 1:num_dir
%         MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
%         RT_S_data = RT_Smean(cd).(strcat('dir', num2str(ndir)))(1,:);
%         RT_S_xy = find(time-abs(RT_S_data)>0,1,'First');
%         MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);   
% 
%         MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
%         RT_K_data = RT_Kmean(cd).(strcat('dir', num2str(ndir)))(1,:);
%         RT_K_xy = find(time-abs(RT_K_data)>0,1,'First');
%         MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);   
%         % RT cycle
%         subplot(2,4,ndir)
%         plot(MS_data_x(1:RT_S_xy), MS_data_y(1:RT_S_xy), 'r'); % RT mean dir
%         hold on;
%         plot(MS_data_x(RT_S_xy+1:end), MS_data_y(RT_S_xy+1:end), 'b');
%         hold on;
%         plot(MK_data_x(1:RT_K_xy), MK_data_y(1:RT_K_xy), 'k'); % RT mean dir
%         hold on;
%         plot(MK_data_x(RT_K_xy+1:end), MK_data_y(RT_K_xy+1:end), 'g');
%         hold on;
%         labeltitle = strcat('Dir ', num2str(ndir));
%         title(labeltitle)
%         label_text = strcat('Dir ', num2str(ndir)); 
%         text(MS_data_x(end), MS_data_y(end), label_text, ...
%             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
%         text(MK_data_x(end), MK_data_y(end), label_text, ...
%             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'g');
%             % Creazione di oggetti fittizi per la legenda
%             h1 = plot(nan, nan, 'b');
%             h2 = plot(nan, nan, 'g');
% 
%             % Aggiunta della legenda
%             legend([h1 h2], {'Monkey - S', 'Monkey - K'}, 'Location', 'best');
%     end
%     sgtitle(append('Condition -',' ',condition_name{cd}),'Color',color_name{cd})
%      XY_RTsubplot_path = strcat(dir_path,'\XY_RTsubplot');
%     if ~exist(XY_RTsubplot_path,'dir')
%         mkdir(XY_RTsubplot_path)
%     end
%     name_fig = fullfile(XY_RTsubplot_path,strcat(condition_name{cd},'XY_RTsubplot.png'));
%     saveas(gcf,name_fig,'png');
% end

%% plot cursor x,y position without RT
if exist('plot_withoutRT','var')
    time = linspace(0,1,1000);
    for cd=1:num_cond
        figure
        for ndir = 1:num_dir
            MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
            RT_S_data = RT_Smean(cd).(strcat('dir', num2str(ndir)))(1,:);
            RT_S_xy = find(time-abs(RT_S_data)>0,1,'First');
            MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);

            MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
            RT_K_data = RT_Kmean(cd).(strcat('dir', num2str(ndir)))(1,:);
            RT_K_xy = find(time-abs(RT_K_data)>0,1,'First');
            MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);
            % RT cycle
            plot(MS_data_x(RT_S_xy+1:end), MS_data_y(RT_S_xy+1:end), 'b');
            hold on;
            plot(MK_data_x(RT_K_xy+1:end), MK_data_y(RT_K_xy+1:end), 'g');
            label_text = strcat('Dir ', num2str(ndir));
            text(MS_data_x(end), MS_data_y(end), label_text, ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
            text(MK_data_x(end), MK_data_y(end), label_text, ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'g');
            % Creazione di oggetti fittizi per la legenda
            h1 = plot(nan, nan, 'b');
            h2 = plot(nan, nan, 'g');

            % Aggiunta della legenda
            legend([h1 h2], {'Monkey - S', 'Monkey - K'}, 'Location', 'best');
        end
        Positionxyname = append('Condition -',' ',condition_name{cd});
        title(Positionxyname,'Color',color_name{cd})
        XY_noRT_path = strcat(dir_path,'\XY_noRT');
        if ~exist(XY_noRT_path,'dir')
            mkdir(XY_noRT_path)
        end
        name_fig = fullfile(XY_noRT_path,strcat(condition_name{cd},'XY_noRT.png'));
        saveas(gcf,name_fig,'png');
    end
end
% %% plot cursor x,y position without RT in subplot 
% time = linspace(0,1,1000);
% for cd=1:num_cond
%     figure
%     for ndir = 1:num_dir
%         MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
%         RT_S_data = RT_Smean(cd).(strcat('dir', num2str(ndir)))(1,:);
%         RT_S_xy = find(time-abs(RT_S_data)>0,1,'First');
%         MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);   
% 
%         MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
%         RT_K_data = RT_Kmean(cd).(strcat('dir', num2str(ndir)))(1,:);
%         RT_K_xy = find(time-abs(RT_K_data)>0,1,'First');
%         MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);   
%         % RT cycle
%         subplot(2,4,ndir)
%         plot(MS_data_x(RT_S_xy+1:end), MS_data_y(RT_S_xy+1:end), 'b');
%         hold on;
%         plot(MK_data_x(RT_K_xy+1:end), MK_data_y(RT_K_xy+1:end), 'g');
%         hold on;
%         labeltitle = strcat('Dir ', num2str(ndir));
%         title(labeltitle)
%         label_text = strcat('Dir ', num2str(ndir)); 
%         text(MS_data_x(end), MS_data_y(end), label_text, ...
%             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
%         text(MK_data_x(end), MK_data_y(end), label_text, ...
%             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'g');
%             % Creazione di oggetti fittizi per la legenda
%             h1 = plot(nan, nan, 'b');
%             h2 = plot(nan, nan, 'g');
% 
%             % Aggiunta della legenda
%             legend([h1 h2], {'Monkey - S', 'Monkey - K'}, 'Location', 'best');
%     end
%     sgtitle(append('Condition -',' ',condition_name{cd}),'Color',color_name{cd})
%      XY_noRTsubplot_path = strcat(dir_path,'\XY_noRTsubplot');
%     if ~exist(XY_noRTsubplot_path,'dir')
%         mkdir(XY_noRTsubplot_path)
%     end
%     name_fig = fullfile(XY_noRTsubplot_path,strcat(condition_name{cd},'XY_noRTsubplot.png'));
%     saveas(gcf,name_fig,'png');
% end

%% Velocity profile

%% gradient application
for cd=1:num_cond
    figure
    for ndir = 1:num_dir
        MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
        MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);

        MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
        MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);

        t = linspace(0,1,1000);
        dt = gradient(t);
        dx_s = gradient(MS_data_x);
        dy_s = gradient(MS_data_y);
        dvx_s = dx_s./dt;
        dvy_s = dy_s./dt;

        dv_s = sqrt(dvx_s.^2 + dvy_s.^2);

        dx_k = gradient(MK_data_x);
        dy_k = gradient(MK_data_y);
        dvx_k = dx_k./dt;
        dvy_k = dy_k./dt;

        dv_k = sqrt(dvx_k.^2 + dvy_k.^2);
        subplot(2,4,ndir)
        plot(t, dv_s, 'b');
        hold on;
        plot(t, dv_k, 'g');
        labeltitle = strcat('Dir ', num2str(ndir));
        title(labeltitle)
        legend({'Velocity profile - S','Velocity profile - K'},'Location','northwest');
        xlabel('time [s]')
        ylabel('Velocity Profile []')
    end
    sgtitle(append('Condition -',' ',condition_name{cd}),'Color',color_name{cd})
     VelocityGradient_path = strcat(dir_path,'\VelocityGradient');
    if ~exist(VelocityGradient_path,'dir')
        mkdir(VelocityGradient_path)
    end
    name_fig = fullfile(VelocityGradient_path,strcat(condition_name{cd},'VelocityGradient.png'));
    saveas(gcf,name_fig,'png');
end

%% diff and Smoothing application
for cd=1:num_cond
    figure  
    for ndir = 1:num_dir
        MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
        MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);

        MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
        MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);

        t = linspace(0,1,1000);
        dt = diff(t);

        MS_dx = diff(MS_data_x);
        MS_dy = diff(MS_data_y);

        MS_vx = MS_dx./dt;
        MS_vy = MS_dy./dt;
        MS_velocity = sqrt(MS_vx.^2 + MS_vy.^2);
    
        windowSize = 50;
        ofilter = 3; % Filter order Savitzky-Golay
        %% moving average filter
        MS_smoothMA = movmean(MS_velocity, windowSize);
        %% Savitzky-Golay filter
        MS_smoothSG = sgolayfilt(MS_velocity, ofilter, windowSize+1);

        MK_dx = diff(MK_data_x);
        MK_dy = diff(MK_data_y);

        MK_vx = MK_dx./dt;
        MK_vy = MK_dy./dt;
        MK_velocity = sqrt(MK_vx.^2 + MK_vy.^2);
        %% moving average filter
        MK_smoothMA = movmean(MK_velocity, windowSize);
        %% Savitzky-Golay filter
        MK_smoothSG = sgolayfilt(MK_velocity, ofilter, windowSize+1);
        
        % figure subplot
        subplot(2,4,ndir)
        plot(t(1:end-1), MS_smoothSG, 'b');
        hold on;
        plot(t(1:end-1), MK_smoothSG, 'g');
        labeltitle = strcat('Dir ', num2str(ndir)); 
        title(labeltitle)
        legend({'Velocity profile - S','Velocity profile - K'},'Location','northwest');
        xlabel('time [s]')
        ylabel('Velocity Profile []')
        % text(MS_smoothSG(end), t(end), label_text, ...
        %     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
        % text(MK_smoothSG(end), t(end), label_text, ...
        %     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'g');
    end
    sgtitle(append('Condition -',' ',condition_name{cd}),'Color',color_name{cd})
     Velocitydiff_smooth_path = strcat(dir_path,'\Velocitydiff_smooth');
    if ~exist(Velocitydiff_smooth_path,'dir')
        mkdir(Velocitydiff_smooth_path)
    end
    name_fig = fullfile(Velocitydiff_smooth_path,strcat(condition_name{cd},'Velocitydiff_smooth.png'));
    saveas(gcf,name_fig,'png');
end
