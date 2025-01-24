%% Extract 24 Trials
clear, close all;
rng(10)
%% check if is on server
[~,nn]=system('hostname'); nn=strtrim(nn);
if strcmp(nn,'rk018610')
    isonserver=true;
else
    isonserver=false;
end
%%
session_name                    = 'SK009';%{'SK001','SK009'};  % session name 
idir                            = 1;        % directions -> 1-8 
% fixed parameteres
nLFP                            = 5;        % max number of LFPs
par.BATTAGLIA_CSD               = BATTAGLIA_CSDParams();
S                               = filesep;
if isonserver
    par.BATTAGLIA_CSD.save_dir      = ['~' S 'SAPIENZA' S 'DCM_SERVER' S session_name S 'dir' num2str(idir) S];
else 
    par.BATTAGLIA_CSD.save_dir      = ['~' S 'TESTS' S 'SAPIENZA' S 'DCM_LOCAL' S session_name S 'dir' num2str(idir) S];
end
par.BATTAGLIA_CSD.save_dir = 'D:\main_scriptDCM\BATTAGLIA_Data\';
if ~isfolder(par.BATTAGLIA_CSD.save_dir)
    mkdir(par.BATTAGLIA_CSD.save_dir);
end
fprintf('Saving in %s\n',par.BATTAGLIA_CSD.save_dir);
par.BATTAGLIA_CSD.session_name  = session_name;
par.BATTAGLIA_CSD.whichmodel    = 7;        % Self Model   

try
    dmode;
catch
    dmode = 7;
end
params=getLFPparamsDonnarumma(dmode);
behav_JOINT  = 'H';
behav_EYE    = 'E';
signal_type  = 'Raw';
load([session_name behav_JOINT '_' signal_type]);
Raw_H        = Trials;

[condition, direction, klabels, success] = getClassInfo(Raw_H);
Raw_H       = Raw_H(success);
condition    = condition(success,:);
direction    = direction(success,:);
klabels      = klabels(success,:);

Data = cell(3,8);
for iCond=1:3
    for iDir = 1:8
        Data{iCond,iDir} = find(condition(:,iCond)==1 & direction(:,iDir)==1);
    end
end

n=1;
data24 = NaN(24,1);
for iSel = 1:3
    for iDir = 1:8
        data_sel = Data{iSel,iDir};
        data24(n,1) = data_sel(randi(length(data_sel)));
        n=n+1;
    end
end

Trials = Raw_H(data24);
filepath = 'D:\main_scriptDCM\';
filename = strcat(filepath,'SK009_24H_Raw.mat');
save(filename, 'Trials');