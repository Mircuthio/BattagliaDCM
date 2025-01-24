%% TEST ANOVA per MISURE RIPETUTE
% utile per dati che analizzano stessi soggetti in condizioni o momenti
% temporali differenti


% Dati di esempio per 5 soggetti in 3 condizioni
data = [
    5.2, 4.8, 5.5;  % Soggetto 1
    6.1, 6.0, 5.9;  % Soggetto 2
    5.8, 5.9, 6.1;  % Soggetto 3
    4.9, 5.1, 5.3;  % Soggetto 4
    6.5, 6.3, 6.4   % Soggetto 5
];

% Colonne per ciascuna condizione
% Cond1 Cond2 Cond3

% Creazione della tabella con i dati
subjectIDs = (1:5)';  % Identificativo dei soggetti
tbl = table(subjectIDs, data(:,1), data(:,2), data(:,3), ...
    'VariableNames', {'Subject', 'Cond1', 'Cond2', 'Cond3'});

% Definizione del modello: specifica le variabili per le condizioni
rm = fitrm(tbl, 'Cond1-Cond3 ~ 1', 'WithinDesign', table([1 2 3]', 'VariableNames', {'Condition'}));

% Eseguire l'ANOVA
results = ranova(rm);
disp(results);
