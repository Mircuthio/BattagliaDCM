
function cc_final = circlecenter(distanza_radiale)

cc = NaN(8,2);

cc(1,:) = [0,0];
% Angoli per i cerchi periferici
angoli = linspace(0, 2*pi, 9);  % 8 cerchi disposti in modo radiale
angoli(end) = [];  % Rimuovi l'ultimo angolo per evitare duplicati
for i = 1:length(angoli)
    % Calcola la posizione del centro di ogni cerchio periferico
    x_centro = distanza_radiale * cos(angoli(i));
    y_centro = distanza_radiale * sin(angoli(i));

cc(i+1,:) = [x_centro,y_centro];
end
cc_final(1,:) = cc(1,:); %center
cc_final(2,:) = cc(4,:); %dir1
cc_final(3,:) = cc(3,:); %dir2 
cc_final(4,:) = cc(2,:); %dir3
cc_final(5,:) = cc(9,:); %dir4
cc_final(6,:) = cc(8,:); %dir5
cc_final(7,:) = cc(7,:); %dir6
cc_final(8,:) = cc(6,:); %dir7
cc_final(9,:) = cc(5,:); %dir8