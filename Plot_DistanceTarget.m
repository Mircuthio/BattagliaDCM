% Supponiamo di avere un array target_positions con [x_target, y_target] per ogni direzione
target_positions = [
    0, 5;  % dir1
    -5, 5; % dir2
    -5, 0; % dir3
    -5, -5; % dir4
    0, -5; % dir5
    5, -5; % dir6
    5, 0;  % dir7
    5, 5   % dir8
];

r_dim = 1.81; % Raggio del cerchio target

for cd = 1:num_cond
    figure;
    for ndir = 1:num_dir
        % Centro del bersaglio per questa direzione
        x_target = target_positions(ndir, 1);
        y_target = target_positions(ndir, 2);
        
        % Estrazione dei dati x e y per MS (Monkey-S)
        MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
        MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);
        
        % Estrazione dei dati x e y per MK (Monkey-K)
        MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
        MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);
        
        % Calcolo della distanza dal bersaglio
        MS_distance_target = sqrt((MS_data_x - x_target).^2 + (MS_data_y - y_target).^2);
        MK_distance_target = sqrt((MK_data_x - x_target).^2 + (MK_data_y - y_target).^2);
        
        % Calcolo del tempo (assumendo che i dati siano campionati uniformemente)
        time = 1:length(MS_distance_target); % Modifica se hai un tempo specifico
        
        % Plot della distanza rispetto al bersaglio nel tempo
        subplot(2, 4, ndir); % Disposizione dei subplot
        plot(time, MS_distance_target, 'b', 'LineWidth', 1.5); % Distanza Monkey-S
        hold on;
        plot(time, MK_distance_target, 'g', 'LineWidth', 1.5); % Distanza Monkey-K
        
        % Linea del raggio target
        yline(r_dim, 'r--', 'LineWidth', 1.5); % Linea del raggio del cerchio target
        
        xlabel('Time (samples)');
        ylabel('Distance from target');
        title(strcat('Dir ', num2str(ndir)));
        xlim([1, length(time)]);
        legend({'Monkey-S', 'Monkey-K', 'Target radius'}, 'Location', 'best');
    end
    
    % Testo aggiuntivo con informazioni della condizione
    textString = sprintf("Chamber: %s\nCondition: %d", data_trials(1).Chamber, cd);
    for idir = 1:8
        dirField = strcat('dir', num2str(idir)); 
        textString = sprintf("%s\nD%d-T%d", textString, M_name(cd).(dirField).D, M_name(cd).(dirField).T);
    end
    text(0.95, 0.3, textString, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontSize', 10, 'Color', 'k');
    
    % Imposta la finestra massimizzata
    set(gcf, 'WindowState', 'maximized');
end
