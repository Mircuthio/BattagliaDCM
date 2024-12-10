%% TEST di FRIEDMAN
% test non parametrico su più gruppi DIPENDENTI-APPAIATI
% Dati di esempio per più gruppi (3 condizioni per ogni soggetto)
data = [
    194, 184, 154;
    177, 152, 177;
    212, 177, 198;
    116, 198, 189;
    217, 189, 203;
    287, 203, 172;
    333, 234, 154;
];

% Esegue il test di Friedman
[p, tbl, stats] = friedman(data, 1);

% Output
fprintf('p-value: %.4f\n', p);
