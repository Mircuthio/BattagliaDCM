% Dati di esempio (da sostituire con i tuoi p-value calcolati)
data1 = randn(30, 1);
data2 = randn(30, 1);
data3 = randn(30, 1);
data4 = randn(30, 1);

groups = {data1, data2, data3, data4};
numGroups = numel(groups);

% Matrice dei p-value
p_values = nan(numGroups);

% Calcolo t-test a due a due
for i = 1:numGroups-1
    for j = i+1:numGroups
        [~, p] = ttest2(groups{i}, groups{j});
        p_values(i, j) = p;
    end
end

% Correzione FDR
p_vector = p_values(~isnan(p_values));
fdr_p = mafdr(p_vector, 'BHFDR', true);
corrected_p = nan(size(p_values));
corrected_p(~isnan(p_values)) = fdr_p;

% Soglia di significatività
alpha = 0.05;

% Creazione della heatmap
figure;
h = heatmap(corrected_p, 'CellLabelFormat', '%.4f');
title('P-value corretti (FDR)');
xlabel('Gruppi');
ylabel('Gruppi');

% Rimozione etichetta "NaN" per celle bianche
h.MissingDataLabel = '';

% Definizione della mappa colori personalizzata (bianco per NaN, blu per non significativi, rosso per significativi)
cmap = [1 1 1; 0 0 1; 1 0 0];
colormap(gca, cmap);

% Imposta limiti della scala per colorare solo valori sotto alpha in rosso
h.ColorLimits = [0, alpha]; % Limita i valori significativi

% Imposta il colore bianco per i NaN
h.MissingDataColor = [1 1 1];

% Impostazione della barra dei colori senza personalizzare direttamente l'oggetto colorbar
h.Colormap = cmap;
