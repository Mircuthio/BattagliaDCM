%% TEST di WILCOXON
% test non parametrico su 2 gruppi per dati DIPENDENTI-APPAIATI
% Dati di esempio per due gruppi dipendenti
data1 = [194, 177, 212, 116, 217, 287, 333];
data2 = [184, 152, 177, 198, 189, 203, 172];

% Esegue il test di Wilcoxon per dati appaiati
[p, h,stats] = signrank(data1, data2);

% Output
fprintf('p-value: %.4f\n', p);
fprintf('Hypothesis test result (h): %d\n', h);
