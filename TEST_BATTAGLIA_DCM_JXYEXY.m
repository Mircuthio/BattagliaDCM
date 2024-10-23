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
    test_dir    ='~/TESTS/SAPIENZA/DCM';
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
num_cond = 3; % condizioni 1-SoloS 2-SoloK 3-Joint S-K
condition_name  = {'SoloS';'SoloK';'Joint S-K'};
color_name = {'b','g','y'};

Condition = cell(num_cond,num_dir);
for i=1:num_cond
    for j=1:num_dir
        Condition{i,j} = find([data_trials.Condition]==i & [data_trials.Direction]==j);
    end
end

M = struct(); 
lenM = NaN(num_cond,num_dir);
for i=1:num_cond
    for j=1:num_dir
        M(i).(strcat('dir',num2str(j))) = data_trials(Condition{i,j});
        lenM(i,j) = length(data_trials(Condition{i,j}));
    end
end

%% Cursor Movement (x,y)
Jfieldname = 'JXYEXY';
Jduration = 1000;
MS = struct();
MK = struct();
for i=1:num_cond
    for j=1:num_dir
        M_app = M(i).(strcat('dir',num2str(j)));
        MS_meanX = NaN(lenM(i,j),Jduration);
        MS_meanY = NaN(lenM(i,j),Jduration);
        MK_meanX = NaN(lenM(i,j),Jduration);
        MK_meanY = NaN(lenM(i,j),Jduration);        
        for k = 1:lenM(i,j)
            MS_meanX(k,:) = M_app(k).(Jfieldname)(1,1:Jduration);
            MS_meanY(k,:) = M_app(k).(Jfieldname)(2,1:Jduration);
            MK_meanX(k,:) = M_app(k).(Jfieldname)(5,1:Jduration);
            MK_meanY(k,:) = M_app(k).(Jfieldname)(6,1:Jduration);
        end
        MS(i).(strcat('dir',num2str(j))) = [mean(MS_meanX,1);mean(MS_meanY,1)];
        MK(i).(strcat('dir',num2str(j))) = [mean(MK_meanX,1);mean(MK_meanY,1)];
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
for i=1:num_cond
    for j=1:num_dir
        M_app = M(i).(strcat('dir',num2str(j)));
        RT_Sdir = NaN(lenM(i,j),1);
        RT_Kdir = NaN(lenM(i,j),1);
        for k = 1:lenM(i,j)
            RT_Sdir(k,1) = M_app(k).(strcat(RTfieldname,'_S'));
            RT_Kdir(k,1) = M_app(k).(strcat(RTfieldname,'_K'));
        end
        RT_S(i).(strcat('dir',num2str(j))) = RT_Sdir;
        RT_K(i).(strcat('dir',num2str(j))) = RT_Kdir;
        RT_Smean(i).(strcat('dir',num2str(j))) = mean(RT_Sdir,'omitmissing');
        RT_Kmean(i).(strcat('dir',num2str(j))) = mean(RT_Kdir,'omitmissing');
        RT_Sstd(i).(strcat('dir',num2str(j))) = std(RT_Sdir,'omitmissing');
        RT_Kstd(i).(strcat('dir',num2str(j))) = std(RT_Kdir,'omitmissing');
    end
end

%% RT barplot
    RT_SplotMean = (cell2mat(struct2cell(RT_Smean')))';
    RT_KplotMean = (cell2mat(struct2cell(RT_Kmean')))';
    RT_SplotStd = (cell2mat(struct2cell(RT_Sstd')))';
    RT_KplotStd = (cell2mat(struct2cell(RT_Kstd')))';
for cd = 1:num_cond
    RT_plotMean = 1000*[RT_SplotMean(cd,:);RT_KplotMean(cd,:)];
    RT_plotStd = 1000*[RT_SplotStd(cd,:);RT_KplotStd(cd,:)];
    categories = {'dir1', 'dir2', 'dir3', 'dir4', 'dir5', 'dir6', 'dir7', 'dir8'};
    figure
    RT_barplot = bar(RT_plotMean');
    set(gca, 'XTickLabel', categories); 
    legend({'Monkey - S', 'Monkey - K'}, 'Location', 'northwest'); 
    set(gca, 'XTickLabelRotation', 45); 
    RT_barplot(1).FaceColor = 'b'; 
    RT_barplot(2).FaceColor = 'g'; 
    RTtitle_name = append('Reaction Time -',' ',condition_name{cd});
    ylabel('RT [ms]');
    title(RTtitle_name,'Color',color_name{cd});
    hold on;
    x = (1:length(categories));
    e1 = errorbar(x - 0.15, RT_plotMean(1,:), RT_SplotStd(1,:), ...
                  'vertical', 'linestyle', 'none', 'CapSize', 10, 'Color', 'black', 'LineWidth', 1.5);
    e1.Annotation.LegendInformation.IconDisplayStyle = 'off'; 

    e2 = errorbar(x + 0.15, RT_plotMean(2,:), RT_SplotStd(2,:), ...
                  'vertical', 'linestyle', 'none', 'CapSize', 10, 'Color', 'black','LineWidth', 1.5);
    e2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    hold off;
end

%% RT boxplot
for cd = 1:num_cond
    fieldsRT_S = fieldnames(RT_S);
    tempRT_S = cell(numel(fieldsRT_S), 1);
    for i = 1:numel(fieldsRT_S)
        tempRT_S{i} = RT_S(cd).(fieldsRT_S{i});
    end
    RT_Sbox = vertcat(tempRT_S{:});

    fieldsRT_K = fieldnames(RT_K);
    tempRT_K = cell(numel(fieldsRT_K), 1);
    for i = 1:numel(fieldsRT_K)
        tempRT_K{i} = RT_K(cd).(fieldsRT_K{i});
    end
    RT_Kbox = vertcat(tempRT_K{:});
    groups = [repmat({'S-dir1'}, lenM(cd, 1), 1); repmat({'S-dir2'}, lenM(cd, 2), 1);
        repmat({'S-dir3'}, lenM(cd, 3), 1); repmat({'S-dir4'}, lenM(cd, 4), 1);
        repmat({'S-dir5'}, lenM(cd, 5), 1); repmat({'S-dir6'},lenM(cd, 6), 1);
        repmat({'S-dir7'}, lenM(cd, 7), 1); repmat({'S-dir8'}, lenM(cd, 8), 1);
        repmat({'K-dir1'}, lenM(cd, 1), 1); repmat({'K-dir2'}, lenM(cd, 2), 1);
        repmat({'K-dir3'}, lenM(cd, 3), 1); repmat({'K-dir4'}, lenM(cd, 4), 1);
        repmat({'K-dir5'}, lenM(cd, 5), 1); repmat({'K-dir6'}, lenM(cd, 6), 1);
        repmat({'K-dir7'}, lenM(cd, 7), 1); repmat({'K-dir8'}, lenM(cd, 8), 1)];

    categories = {'dir1', 'dir2', 'dir3', 'dir4', 'dir5', 'dir6', 'dir7', 'dir8'};
    figure
    boxplot([RT_Sbox;RT_Kbox], groups, 'Colors', [0 0 1; 1 1 0]);
    set(gca, 'XTickLabel', categories);
    legend({'Monkey - S', 'Monkey - K'}, 'Location', 'northwest');
    set(gca, 'XTickLabelRotation', 45);
    RTtitle_name = append('Reaction Time -',' ',condition_name{cd});
    ylabel('RT [ms]');
    title(RTtitle_name,'Color',color_name{cd});
end


    %% plot x,y position
for cd=1:num_cond
    figure
    for dir = 1:num_dir
        MS_data_x = MS(cd).(strcat('dir', num2str(dir)))(1,:);
        MS_data_y = MS(cd).(strcat('dir', num2str(dir)))(2,:);

        MK_data_x = MK(cd).(strcat('dir', num2str(dir)))(1,:);
        MK_data_y = MK(cd).(strcat('dir', num2str(dir)))(2,:);

        plot(MS_data_x, MS_data_y, 'b');
        hold on;
        plot(MK_data_x, MK_data_y, 'g');
        label_text = strcat('Dir ', num2str(dir)); 
        text(MS_data_x(end), MS_data_y(end), label_text, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
        text(MK_data_x(end), MK_data_y(end), label_text, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'g');
    end
    title(strcat('Condition ',condition_name{cd}),'Color',color_name{cd})
end

        %% Velocity profile
for cd=1:num_cond
    figure
    for dir = 1:num_dir
        MS_data_x = MS(cd).(strcat('dir', num2str(dir)))(1,:);
        MS_data_y = MS(cd).(strcat('dir', num2str(dir)))(2,:);

        MK_data_x = MK(cd).(strcat('dir', num2str(dir)))(1,:);
        MK_data_y = MK(cd).(strcat('dir', num2str(dir)))(2,:);

        t = linspace(0,1,1000);
        dt = diff(t);

        MS_dx = diff(MS_data_x);
        MS_dy = diff(MS_data_y);
        tdiff = diff(t);
        MS_vx = MS_dx./tdiff;
        MS_vy = MS_dy./tdiff;
        MS_velocity = sqrt(MS_vx.^2 + MS_vy.^2);
    
        windowSize = 50;
        ofilter = 3; % Filter order Savitzky-Golay
        %% moving average filter
        MS_smoothMA = movmean(MS_velocity, windowSize);
        %% Savitzky-Golay filter
        MS_smoothSG = sgolayfilt(MS_velocity, ofilter, windowSize+1);

        MK_vx = MK_dx./tdiff;
        MK_vy = MK_dy./tdiff;
        MK_velocity = sqrt(MK_vx.^2 + MK_vy.^2);
        %% moving average filter
        MK_smoothMA = movmean(MK_velocity, windowSize);
        %% Savitzky-Golay filter
        MK_smoothSG = sgolayfilt(MK_velocity, ofilter, windowSize+1);


    plot(t(1:end-1), MS_velocity, 'b');
    hold on;
    plot(t(1:end-1), MK_velocity, 'g');
        label_text = strcat('Dir ', num2str(dir)); 
        text(MS_velocity(end), t(end), label_text, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'black');
        text(MK_velocity(end), t(end), label_text, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'red');
    end
    title(strcat('Condition ',condition_name{cd}),'Color',color_name{cd})
end

%% Step 1: CSD compute
fprintf('Step 1: csd compute\n');
signal_process                          = 'CSD';
par.csdCompute                          = csdComputeParams();            
par.csdCompute.InField                  = 'LFP';
par.csdCompute.OutField                 = signal_process;
par.csdCompute.rsfactor                 = 0.2;
par.csdCompute.optrescale               = 3;%-1;
data_trials                             = csdCompute(data_trials,par.csdCompute);
test_dir                                = [test_dir '_opt' num2str(par.csdCompute.optrescale) S];





%% learning
fprintf('Step 2: Monkey S-K model\n')
par.dcmJointModel.whichmodel            = 5;
par.dcmJointModel.custom_model          = [];%'customPriorsv3';
par.dcmJointModel.InField               = signal_process;
par.dcmJointModel.donlfp                = false;
par.dcmJointModel.isdemo                = getSelectionIndexes(idir,1:3); % get S, K and S-K Trials
par.dcmJointModel.S_PRIORS              = [];%DCM_S{S_winner}.fn;
par.dcmJointModel.K_PRIORS              = [];%DCM_K{K_winner}.fn;
out.dcmJointModel.DCM                   = dcmJointModel(data_trials,par.dcmJointModel);





%% save test 
% custom directory and file name string construction
[~,Labels]  = getJointMonkeysLabels(1:24);
iConds       = find(ismember(Labels,DCM_Joint.xU.name));
save_dir    =[test_dir session_name S];
save_dir    =[save_dir 'Dir'];
for ind=1:length(idir)
    save_dir=[save_dir '_' num2str(idir(ind))];
end
save_dir    =[save_dir S];
save_file   =['LFP_M' num2str(par.dcmJointModel.whichmodel) '_' session_name];
save_file   =[save_file '_D' num2str(par.BattagliaArrangeTrials.dmode)];
save_file   =[save_file '_S' num2str(par.BattagliaArrangeTrials.selS)];
save_file   =[save_file '_K' num2str(par.BattagliaArrangeTrials.selK)];

save_file=[save_file '_C'];
for inc=1:length(iConds)
    save_file=[save_file '_' num2str(iConds(inc))];
end
if ~isempty (par.dcmJointModel.custom_model)
    save_file=[save_file '_' par.dcmJointModel.custom_model];
end
if par.dcmJointModel.donlfp
    save_file = [save_file '_DONLFP'];
end
fprintf('Saving in %s\n',[save_dir save_file]);
DCM         =out.dcmJointModel.DCM;
DCM.file    =save_file;
DCM.fullpath=[save_dir, save_file];
DCM.par     =par.dcmJointModel;
if ~isfolder(save_dir)
    mkdir(save_dir);
end
save(DCM.fullpath,'DCM','data_trials');
%% stop if is on server
if isonserver; return; end
%% LFP
iTrial          = 1;
fprintf('Showing Trial %g\n',iTrial);
hfg.lfp         = figure;
hold on; box on; grid on;
lw              = 2;
InField         = 'LFP';
xfld            = 't';
TField          = [xfld InField];
col1            = [0,0,0];
nSources    = size(data_trials(iTrial).(InField),1);
tiledlayout(nSources,1);
sigtime     = data_trials(iTrial).(TField);
sigtime     = sigtime-sigtime(1);
for iSource=1:nSources
    nexttile;
    hold on; box on; grid on;
    plot(sigtime,data_trials(iTrial).(InField)(iSource,:),'color',col1,'linewidth',lw);
    xlabel('time [s]')
    ylabel('LFP [mV]') 
    title(['Source ' num2str(iSource)]);
end
sgtitle(['Low Field Potential, Trial ' num2str(iTrial)])
%% CSD
iTrial          = 2;
fprintf('Showing CSD result of Trial %g\n',iTrial);
hfg.cross       = figure;
col1            = [1,0,0];
col2            = [0,1,0];
lw              = 3;
InField         = 'CSD';
xfld            = 't';
TField          = [xfld,InField];
[~,nRows,nCols] = size(data_trials(iTrial).(InField));
p               = [0,0,1500,600];
fs              = 13;
tiledlayout(nRows,nCols);
for iRow=1:nRows
    for iCol=1:nCols
        nexttile;
        hold on; box on; grid on;
        % plot original
        xl  = [data_trials(iTrial).(TField)(1),data_trials(iTrial).(TField)(end)];
        plot(data_trials(iTrial).(TField), ...
             real(data_trials(iTrial).(InField)(:,iRow,iCol)),'color',col1,'linewidth',lw,'linestyle',':');%,M.Hz,real(CSD(:,1,2)),':')
        % plot reconstructed
        plot(data_trials(iTrial).(TField), ...
             real(DCM.Hc{iTrial}(:,iRow,iCol)),'color',col2,'linewidth',lw);%,M.Hz,real(CSD(:,1,2)),':')
        xlim(xl)
        title(['Source ' num2str(iRow) ' vs Source' num2str(iCol)])
        if iCol==nCols && iRow==nRows
            legend('original','reconstructed');
        end
        if iRow==nRows
            xlabel('Freq [Hz]');
        end
        if iCol==1
            ylabel('CSD [dB/Hz]');
        end
        set(gca,'fontsize',fs);
    end
end
trialError  = computeError(data_trials(iTrial).(InField),DCM.Hc{iTrial});
st          = sprintf('[cross]-spectral density - RMSE: %0.3f', trialError.RMSE);
sgtitle(st,'fontsize',fs+2)
set(hfg.cross,'Position',p);
set(hfg.cross, 'PaperPositionMode','auto');

%%
% plot parameters and estimates
%--------------------------------------------------------------------------
hfg.expects = figure;
hold on; box on; grid on;
bar(exp(spm_vec(DCM.Ep)))
title('conditional expectation')
errorbar(exp(spm_vec(DCM.Ep)),exp(DCM.Cp))

%% some plots
plot_spm_dcm_csd_results(DCM,'Coupling (B)',figure);
plot_spm_dcm_csd_results(DCM,'trial-specific effects',figure);

% DCM_RESULTS_LFP_KS(DCM);
%%
