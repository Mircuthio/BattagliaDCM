%% Omogenizza i dati

function data_out = omogdata(data)

len = NaN(length(data),1);
for ns=1:length(data)
    len(ns) = length(data{ns});
end
lenmin = min(len);
data_out = struct();
for ns=1:length(data)
    data_app = data{ns};
    if length(data_app) > lenmin
        ndiff = length(data_app) - lenmin;
        idx = randperm(length(data_app),ndiff);
        data_app(idx) = [];
    end
    data_out(ns).group = data_app;
end
data_out = cell2mat(struct2cell(data_out'));