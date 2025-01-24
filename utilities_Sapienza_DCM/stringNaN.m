% create stringNaN

function infoNaN = stringNaN(data,num_cond,num_dir)
latexString = struct();
for i=1:num_cond
    for j=1:num_dir
        data_app = data(i).(strcat('dir',num2str(j)));
        latexString(i).(strcat('dir',num2str(j))) = ['_{' num2str(data_app.S) '}' '{' num2str(data_app.len) '_' '{' num2str(data_app.K) '}}'];
        % latexString(i).(strcat('dir',num2str(j))) = [num2str(data_app.S) '_{' num2str(data_app.len) '}^{' num2str(data_app.K) '}$'];

    end
end
infoNaN = struct2cell(latexString')';
