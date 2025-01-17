for cd = 1:num_cond
    figure;
    for ndir = 1:num_dir
        % Estrazione dei dati x e y per MS (Monkey-S)
        MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
        MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);
        
        % Estrazione dei dati x e y per MK (Monkey-K)
        MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
        MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);
        
        % Calcolo della lunghezza della traiettoria (Monkey-S)
        dx_MS = diff(MS_data_x);
        dy_MS = diff(MS_data_y);
        traj_length_MS = sum(sqrt(dx_MS.^2 + dy_MS.^2));

        % Calcolo della lunghezza della traiettoria (Monkey-K)
        dx_MK = diff(MK_data_x);
        dy_MK = diff(MK_data_y);
        traj_length_MK = sum(sqrt(dx_MK.^2 + dy_MK.^2));

        % Calcolo della distanza rispetto al centro (0,0) per MS e MK
        MS_distance = sqrt(MS_data_x.^2 + MS_data_y.^2);
        MK_distance = sqrt(MK_data_x.^2 + MK_data_y.^2);
        
        % Calcolo del tempo (assumendo che i dati siano campionati uniformemente)
        time = 1:length(MS_distance); % Modifica se hai un tempo specifico
        ticktime = 1000*Time(cd).(strcat('dir',num2str(ndir)));
        % Plot delle distanze nel tempo
        subplot(2, 4, ndir); % Disposizione dei subplot
        plot(time, MS_distance, 'b', 'LineWidth', 1.5); % Plot distanza per Monkey-S
        hold on;
        plot(time, MK_distance, 'g', 'LineWidth', 1.5); % Plot distanza per Monkey-K
        xlabel('Time (ms)');
        ylabel('Distance from center');
        title(strcat('Dir ', num2str(ndir)));
        xlim([1, length(time)]);
        timestep = round(floor(length(ticktime)/15));
        xticks(round(ticktime(1:timestep:end)))
        legend({'Monkey-S', 'Monkey-K'}, 'Location', 'best');
    
        text(8, 6, sprintf('Monkey-S = %.1f \n Monkey-K = %.1f', ...
            traj_length_MS, traj_length_MK), ...
            'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'left');
    % % Testo aggiuntivo con le informazioni della condizione e dei parametri
    % textString = sprintf("Chamber: %s\nCondition: %d", data_trials(1).Chamber, cd);
    % for idir = 1:8
    %     dirField = strcat('dir', num2str(idir)); 
    %     textString = sprintf("%s\nD%d-T%d", textString, M_name(cd).(dirField).D, M_name(cd).(dirField).T);
    % end
    % text(0.95, 0.3, textString, 'Units', 'normalized', ...
    %     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    %     'FontSize', 10, 'Color', 'k');
    end
    sgtitle(condition_name{cd}, 'Color', color_name{cd});
    set(gcf, 'WindowState', 'maximized');
end
