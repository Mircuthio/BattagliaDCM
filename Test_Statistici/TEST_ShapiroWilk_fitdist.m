%% TEST di NORMALITA tramite fitdist analogo matlab di Test di Shapiro-Wilk (swtest):

% Gruppo di dati
data = [10, 12, 15, 13, 14];

% occorre eliminare i NaN
data_clean = data(~isnan(data));

[h, p] = swtest(data_clean);


% Adatta una distribuzione normale ai dati
pd = fitdist(data', 'Normal');

% Mostra i parametri della distribuzione adattata
disp(['Media: ', num2str(pd.mu)]);
disp(['Deviazione standard: ', num2str(pd.sigma)]);
% Gruppo di dati
data = [10, 12, 15, 13, 14];

% Adatta la distribuzione normale
pd = fitdist(data', 'Normal');

% Grafico Q-Q per la distribuzione normale
qqplot(data, pd);
% Test di Kolmogorov-Smirnov per verificare la normalità
[h, p] = kstest((data - pd.mu) / pd.sigma); % Normalizza i dati
disp(['p-value del test di normalità: ', num2str(p)]);
