% Dati per tre gruppi
group1 = [10, 12, 15, 13, 14];
group2 = [22, 25, 28, 26, 27];
group3 = [18, 20, 23, 21, 24];

% Concatenare i dati in un'unica matrice
data = [group1, group2, group3];

%% occorre eliminare i NaN
data_clean = data(~any(isnan(data), 2), :);

% Creare un vettore di gruppi
group = [ones(1, length(group1)), 2*ones(1, length(group2)), 3*ones(1, length(group3))];

% Test di Kruskal-Wallis
[p, tbl, stats] = kruskalwallis(data, group);

% Visualizzare i risultati
disp(p);
disp(tbl);
disp(stats);
