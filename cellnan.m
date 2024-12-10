%% NaN in cell

function data_clean = cellnan(data)

data_clean = cell(size(data));
for i=1:length(data)
    data_app = data{i};
    data_clean{i} = data_app(~isnan(data_app));
end
