%% TEst di levene
% Gruppi di dati
group1 = [10, 12, 15, 13, 14];
group2 = [20, 22, 18, 21, 19, 23];
group3 = [30, 32, 31, 28, 35];

%% occorre eliminare i NaN
data_clean = data(~any(isnan(data), 2), :);  % Rimuove le righe con NaN

% Esegui il test di Levene
p_value = Levenetest({group1, group2, group3});

% Mostra il p-value
disp(['p-value del test di Levene: ', num2str(p_value)]);
