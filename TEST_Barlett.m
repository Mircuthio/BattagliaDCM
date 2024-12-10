% Gruppi di dati
group1 = [10, 12, 15, 13, 14];
group2 = [20, 22, 18, 21, 19, 23];
group3 = [30, 32, 31, 28, 35];

% Ristruttura i dati in una matrice (ogni colonna è un gruppo)
data = [group1', group2', group3'];  % Trasposta per garantire che ogni gruppo sia una colonna
%% occorre eliminare i NaN
data_clean = data(~any(isnan(data), 2), :);  % Rimuove le righe con NaN

% Esegui il test di Bartlett
[p_value, stats] = bartlett(data);

% Mostra il p-value
disp(['p-value del test di Bartlett: ', num2str(p_value)]);
disp('Statistiche del test di Bartlett:');
disp(stats);

