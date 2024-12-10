% function addcircle(raggio_cerchi)

function  addcircle(raggio_cerchi)

distanza_radiale = 8;  % Distanza radiale per i cerchi periferici
theta = linspace(0, 2*pi, 100);  % Angoli per il contorno dei cerchi

% Cerchio centrale (contorno rosso)
x_centrale = raggio_cerchi * cos(theta);
y_centrale = raggio_cerchi * sin(theta);
plot(x_centrale, y_centrale, 'r', 'LineStyle','--','LineWidth', 0.5);  % Contorno rosso, senza riempimento

% Angoli per i cerchi periferici
angoli = linspace(0, 2*pi, 9);  % 8 cerchi disposti in modo radiale
angoli(end) = [];  % Rimuovi l'ultimo angolo per evitare duplicati

% Disegna i cerchi periferici (contorni neri)
for i = 1:length(angoli)
    % Calcola la posizione del centro di ogni cerchio periferico
    x_centro = distanza_radiale * cos(angoli(i));
    y_centro = distanza_radiale * sin(angoli(i));
    
    % Calcola i punti per il contorno del cerchio periferico
    x_periferico = x_centro + raggio_cerchi * cos(theta);
    y_periferico = y_centro + raggio_cerchi * sin(theta);
    
    % Disegna il contorno nero per ogni cerchio periferico
    plot(x_periferico, y_periferico, 'k','LineStyle','--', 'LineWidth', 0.5);  % Contorno nero
end
