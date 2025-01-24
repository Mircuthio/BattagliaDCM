% Grafico distanza da un target

function plotTargetDistance(data1,data2,Time,M_name,par)

condition_name = par.condition_name;
color_name = par.color_name;
chamber = par.chamber;


target_positions = par.circlec;
target_position_names = ['C','D1','D2','D3','D4','D5','D6','D7','D8'];
MS = data1;
MK = data2;
num_cond = length(MS);
num_dir = length(fieldnames(MS));
for cd = 1:num_cond
    for ndir = 1:num_dir
        figure;
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
        target = ['Dir ', num2str(ndir)];
        for ntarg = 1:length(target_positions)
            x_target = target_positions(ntarg, 1);
            y_target = target_positions(ntarg, 2);
            if strcmp(target_position_names(ntarg),'C')
                % Calcolo della distanza rispetto al centro (0,0) per MS e MK
                MS_distance = sqrt(MS_data_x.^2 + MS_data_y.^2);
                MK_distance = sqrt(MK_data_x.^2 + MK_data_y.^2);
                y_max = max(max(MS_distance),max(MK_distance));
                traj_x = 12;
                traj_y = 6;
                text_x = 0.95;
                text_y = 0.95;
                y_labelname = strcat('|| Jxy (t) - C_{0}||');
                title_name = strcat(condition_name{cd},' \newline Target Dir',num2str(ndir));
                name_fig = strcat(condition_name{cd},'_XY_Dist_',target,'Center');
            else
                MS_distance = sqrt((MS_data_x - x_target).^2 + (MS_data_y - y_target).^2);
                MK_distance = sqrt((MK_data_x - x_target).^2 + (MK_data_y - y_target).^2);
                y_max = max(max(MS_distance),max(MK_distance));
                traj_x = 12;
                traj_y = 1;
                text_x = 0.15;
                text_y = 0.15;
                y_labelname = strcat('|| Jxy (t) - C_{Dir',num2str(ntarg-1),'}||');
                title_name = strcat(condition_name{cd},' \newline Target Dir',num2str(ndir));
                name_fig = strcat(condition_name{cd},'_XY_Dist_',target,'_Dir_',num2str(ntarg-1));
            end

            % Calcolo del tempo (assumendo che i dati siano campionati uniformemente)
            time = 1:length(MS_distance); % Modifica se hai un tempo specifico
            ticktime = 1000*Time(cd).(strcat('dir',num2str(ndir)));
            % Plot delle distanze nel tempo
            plot(time, MS_distance, 'b', 'LineWidth', 1.5); % Plot distanza per Monkey-S
            hold on;
            plot(time, MK_distance, 'g', 'LineWidth', 1.5); % Plot distanza per Monkey-K
            xlabel('Time (ms)');
            ylabel(y_labelname);
            title(title_name, 'Interpreter','tex','Color', color_name{cd})
            xlim([1, length(time)])
            ylim([0 y_max+y_max/6])
            timextick = 0:50:round(ticktime(end));
            xticks(timextick)
            % timestep = round(floor(length(ticktime)/15));
            % xticks(round(ticktime(1:timestep:end)))
            if par.trjleng == 1
                text(traj_x, traj_y, strcat('Trajectory Length: \newline',sprintf('Monkey-S = %.1f \nMonkey-K = %.1f', ...
                    traj_length_MS, traj_length_MK)), ...
                    'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left');
            end
            % Testo aggiuntivo con le informazioni della condizione e dei parametri
            textString = sprintf("Chamber: %s\nCondition: %d", chamber, cd);
            dirField = strcat('dir', num2str(ndir));
            textString = sprintf("%s\nD%d-T%d", textString, M_name(cd).(dirField).D, M_name(cd).(dirField).T);

            text(text_x, text_y, textString,  'Units', 'normalized',...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                'FontSize', 12, 'Color', 'k');
            legend({'Monkey-S', 'Monkey-K'}, 'Location', 'northwestoutside');
            ax = gca;
            set(ax, 'FontSize', 12);
            set(gcf, 'WindowState', 'maximized');

            XY_path = strcat(par.dir_path,'XY_DistanceFromONE\',target);

            par.savePlotEpsPdfMat.dir_png = strcat(XY_path,'\PNGs\');
            par.savePlotEpsPdfMat.dir_pdf = strcat(XY_path,'\PDFs\');
            par.savePlotEpsPdfMat.dir_mat = strcat(XY_path,'\MATfiles\');

            par.savePlotEpsPdfMat.file_name = name_fig;
            savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
            close all
        end
    end
end
