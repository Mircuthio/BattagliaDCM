function plotJxyCircle(data1,data2,M_name,par)

condition_name = par.condition_name;
color_name = par.color_name;
chamber = par.chamber;
session = par.session;
MS = data1;
MK = data2;
num_cond = length(MS);
num_dir = length(fieldnames(MS));
for cd=1:num_cond
    for ndir = 1:num_dir
        figure
        MS_data_x = MS(cd).(strcat('dir', num2str(ndir)))(1,:);
        MS_data_y = MS(cd).(strcat('dir', num2str(ndir)))(2,:);

        MK_data_x = MK(cd).(strcat('dir', num2str(ndir)))(1,:);
        MK_data_y = MK(cd).(strcat('dir', num2str(ndir)))(2,:);

        plot(MS_data_x, MS_data_y, 'b');
        xlim([-10 10])
        ylim([-10 10])
        hold on;
        plot(MK_data_x, MK_data_y, 'g');
        xlim([-10 10])
        ylim([-10 10])
        hold on;
        label_text = strcat('Dir ', num2str(ndir));
        if ~isnan(MS_data_x)
            text(MS_data_x(end), MS_data_y(end), label_text, ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
        end
        if ~isnan(MK_data_x)
            text(MK_data_x(end), MK_data_y(end), label_text, ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', [0.3, 0.6, 0.4]);
        end
        % Creazione di oggetti fittizi per la legenda
        h1 = plot(nan, nan, 'b');
        h2 = plot(nan, nan, 'g');

        % Aggiungo i cerchi
        r_dim = par.r_dim;
        %% function to add circles
        addcircle(r_dim);
        % Aggiunta della legenda
        legend([h1 h2], {'Monkey - S', 'Monkey - K'}, 'Location', 'northwest');
        % Positionxyname = append('Condition -',' ',condition_name{cd});
        Positionxyname = append(condition_name{cd});
        title(Positionxyname,'Color',color_name{cd})
        textString = sprintf("Session: %s\nChamber: %s\nCondition: %d", session,chamber, cd);
        dirField = strcat('dir', num2str(ndir));
        textString = sprintf("%s\nD%d-T%d", textString, M_name(cd).(dirField).D, M_name(cd).(dirField).T);

        text(0.95, 0.3, textString, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'FontSize', 10, 'Color', 'k');
        set(gcf, 'WindowState', 'maximized');

        XY_path = strcat(par.dir_path,'\XY_CircleCurve');

        par.savePlotEpsPdfMat.dir_png = strcat(XY_path,'\PNGs\');
        par.savePlotEpsPdfMat.dir_pdf = strcat(XY_path,'\PDFs\');
        par.savePlotEpsPdfMat.dir_mat = strcat(XY_path,'\MATfiles\');

        name_fig = strcat(condition_name{cd},'_XY_Dir_',num2str(ndir));
        par.savePlotEpsPdfMat.file_name = name_fig;
        savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
        close all
    end
end