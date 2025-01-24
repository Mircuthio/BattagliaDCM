% function data_final = ArrangeTrialsConcatenate(params)

function data_final = ArrangeTrialsConcatenate(isdemo,session_name)

params.BattagliaArrangeTrials              = BattagliaArrangeTrialsParams();
% par.BattagliaArrangeTrials.whichmodel   = 7;        % 7 Model for action session. 5 Model for all sessions
params.BattagliaArrangeTrials.isdemo       = isdemo;        % getSelectionIndexes(1,1:3);
params.BattagliaArrangeTrials.session_name = session_name;  % which session
ind_NaN = cell(5,5);
data_trials = struct();
for i=1:5
    params.BattagliaArrangeTrials.selS         = i;
    for j=1:5
        params.BattagliaArrangeTrials.selK         = j;
        try
            data_trials(i).(sprintf('K_%d',j)) = BattagliaArrangeTrials(params.BattagliaArrangeTrials);
        catch
            ind_NaN{i,j} = [i,j];
        end
    end
end

data_final = []; 
for i = 1:numel(data_trials)
    fields = fieldnames(data_trials(i)); 
    for j = 1:numel(fields)
        currentField = data_trials(i).(fields{j}); 
        if isempty(data_final)

            data_final = currentField(:);
        else
            data_final = cat(1, data_final, currentField(:));
        end
    end
end
