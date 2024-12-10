%% Grafico Q-Q Quantile-Quantile
% Gruppo di dati
group1 = [10, 12, 15, 13, 14];
group2 = [20, 22, 18, 21, 19, 23];

% Grafico Q-Q per il gruppo1
figure;
subplot(1, 2, 1);
qqplot(group1);
title('Q-Q Plot - Group 1');

% Grafico Q-Q per il gruppo2
subplot(1, 2, 2);
qqplot(group2);
title('Q-Q Plot - Group 2');
