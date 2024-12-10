% Dati di esempio (sostituisci con i tuoi dati)
data = [randn(30, 1); randn(30, 1) + 1; randn(30, 1) + 2; randn(30, 1) + 3];
gruppi = [repmat({'G1'}, 30, 1); repmat({'G2'}, 30, 1); repmat({'G3'}, 30, 1); repmat({'G4'}, 30, 1)];

% Esegui l'ANOVA
[p, tbl, stats] = anova1(data, gruppi, 'off');

% Correzione di Tukey
result = multcompare(stats, 'CType', 'tukey-kramer');
% Numero di gruppi
numGroups = length(stats.gnames);

% Inizializza una matrice vuota per i p-value
p_values_corrected = nan(numGroups);

% Riempie la matrice con i p-value di Tukey
for i = 1:size(result, 1)
    g1 = result(i, 1);  % Indice del gruppo 1
    g2 = result(i, 2);  % Indice del gruppo 2
    p_val = result(i, 6);  % p-value corretto
    p_values_corrected(g1, g2) = p_val;
end
% Visualizza il grafico `multcompare` per i confronti a coppie
figure;
multcompare(stats, 'CType', 'tukey-kramer');
title('Confronti multipli (Tukey)');
% Soglia di significatività
alpha = 0.05;

% Creazione della heatmap
figure;
h = heatmap(p_values_corrected, 'CellLabelFormat', '%.4f');
title('P-value corretti (Tukey)');
xlabel('Gruppi');
ylabel('Gruppi');

% Rimozione etichetta "NaN" per celle bianche
h.MissingDataLabel = '';

% Definizione della mappa colori personalizzata (bianco per NaN, blu per non significativi, rosso per significativi)
cmap = [1 1 1; 0 0 1; 1 0 0];
colormap(gca, cmap);

% Imposta limiti della scala per colorare solo valori sotto alpha in rosso
h.ColorLimits = [0, alpha];

% Imposta il colore bianco per i NaN
h.MissingDataColor = [1 1 1];
