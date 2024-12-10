%% TEST di Mann-Whitney

% Gruppi di dati
group1 = [10, 12, 15, 13, 14];
group2 = [22, 25, 28, 26, 27];

%% occorre eliminare i NaN
data_clean1 = data1(~isnan(data1));
data_clean2 = data2(~isnan(data2));

% Test di Mann-Whitney
[p, h, stats] = ranksum(group1, group2);

% p: p-value, h: risultato del test (0 = no differenze significative, 1 = differenze significative)
