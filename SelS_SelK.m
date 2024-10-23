%% function TEST_BATTAGLIA_DCM prova selS e selK
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
ind = cell(5,5);
for i=1:5
    par.BattagliaArrangeTrials.selS         = i;
    for j=1:5
        par.BattagliaArrangeTrials.selK         = j;
        par.BattagliaArrangeTrials.session_name = session_name;  % which session
        try
            data_trials                             = BattagliaArrangeTrials(par.BattagliaArrangeTrials);
        catch
            ind{i,j} = [i,j];
        end
    end
end