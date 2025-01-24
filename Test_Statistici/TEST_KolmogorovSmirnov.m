%% TEST KOLMOGOROV SMIRNOV
% Gruppo di dati
group1 = [10, 12, 15, 13, 14];
group2 = [20, 22, 18, 21, 19, 23];

%% occorre eliminare i nan e non è sensibile a dati diversi
data_clean = data(~isnan(data));
[h, p] = kstest(data_clean);


% Test di Kolmogorov-Smirnov per la normalità (verifica se i dati seguono una normale)
[h1, p1] = kstest((group1 - mean(group1)) / std(group1)); % Normalizza i dati
[h2, p2] = kstest((group2 - mean(group2)) / std(group2)); % Normalizza i dati

% Risultati
disp(['p-value del test di Kolmogorov-Smirnov per group1: ', num2str(p1)]);
disp(['p-value del test di Kolmogorov-Smirnov per group2: ', num2str(p2)]);
