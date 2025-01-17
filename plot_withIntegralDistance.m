% Distanza ideale (raggio del cerchio periferico)
ideal_distance = 10; % Supponiamo un raggio di 10 unità (modifica secondo il caso reale)
for cd = 1:num_cond
    figure
    for ndir = 1:num_dir
        % Estrazione delle traiettorie
        MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1, :);
        MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2, :);

        MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1, :);
        MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2, :);

        % Calcolo della lunghezza della traiettoria (Monkey-S)
        dx_MS = diff(MS_data_x);
        dy_MS = diff(MS_data_y);
        traj_length_MS = sum(sqrt(dx_MS.^2 + dy_MS.^2));

        % Calcolo della lunghezza della traiettoria (Monkey-K)
        dx_MK = diff(MK_data_x);
        dy_MK = diff(MK_data_y);
        traj_length_MK = sum(sqrt(dx_MK.^2 + dy_MK.^2));

        % Plot delle traiettorie
        plot(MS_data_x, MS_data_y, 'b'); % Traiettoria Monkey-S
        hold on;
        plot(MK_data_x, MK_data_y, 'g'); % Traiettoria Monkey-K
        hold on;

        % Configurazione degli assi
        xlim([-10 10]);
        ylim([-10 10]);

        % Etichette dei punti finali
        label_text = strcat('Dir ', num2str(ndir));
        if ~isnan(MS_data_x)
            text(MS_data_x(end), MS_data_y(end), label_text, ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
                'FontSize', 10, 'Color', 'b');
        end
        if ~isnan(MK_data_x)
            text(MK_data_x(end), MK_data_y(end), label_text, ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
                'FontSize', 10, 'Color', [0.3, 0.6, 0.4]);
        end

        % Aggiunta dei cerchi
        r_dim = 1.81;
        addcircle(r_dim); % Funzione personalizzata

        % Creazione di oggetti fittizi per la legenda
        h1 = plot(nan, nan, 'b');
        h2 = plot(nan, nan, 'g');

        % Aggiunta della legenda
        legend([h1 h2], {'Monkey - S', 'Monkey - K'}, 'Location', 'northwest');

        % Titolo con nome della condizione
        title(condition_name{cd}, 'Color', color_name{cd});

        % Aggiunta della lunghezza della traiettoria nel grafico
        text(-9, 9 - ndir, sprintf('Dir %d: Monkey-S = %.2f, Monkey-K = %.2f', ...
            ndir, traj_length_MS, traj_length_MK), ...
            'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'left');
    end

    % Informazioni aggiuntive sul grafico
    textString = sprintf("Chamber: %s\nCondition: %d", data_trials(1).Chamber, cd);
    for idir = 1:8
        dirField = strcat('dir', num2str(idir));
        textString = sprintf("%s\nD%d-T%d", textString, M_name(cd).(dirField).D, M_name(cd).(dirField).T);
    end
    text(0.95, 0.3, textString, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontSize', 10, 'Color', 'k');

    % Massimizzazione della finestra
    set(gcf, 'WindowState', 'maximized');
end
