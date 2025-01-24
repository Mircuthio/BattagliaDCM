
function RTbarplotSingleDir(RT_plotMean,RT_plotStd,par)

for ndir=1:par.num_dir
    condition_name1 = {'       Monkey S',...
        '       Monkey K',...
        '       Monkey S',...
        '       Monkey K'};
    condition_name2 = {'Act\_S - Obs\_K',...
        'Obs\_S - Act\_K',...
        'Act\_S - Act\_K',...
        'Act\_S - Act\_K'};
    condition_name3 = {'       (Solo S)',...
        '       (Solo K)',...
        '       (Joint S)',...
        '       (Joint K)'};


    % orange triple
    o = {[1 0.5 0],[1 0.7 0.5],[0.8 0.3 0],[1.0, 0.6, 0.0]};
    y = {[1 1 0],[1 1 0.2],[0.8 0.8 0],[1.0, 0.85, 0.0]};
    colorplot = {[0 0 1],[0 1 0],o{2},o{1}};
    colorerror = {[0 1 1],[0.3, 0.8, 0.6],o{4},o{3}};

    figure
    RT_pldirMean = RT_plotMean(:,ndir);
    RT_pldirSt = RT_plotStd(:,ndir);
    RT_barplot = bar(RT_pldirMean,'FaceColor', 'flat');
    for k=1:length(RT_pldirMean)
        RT_barplot.CData(k,:) = colorplot{k};
    end
    % infoNaN_plot = infoNaN(cd,:);
    cat_plot = cellfun(@(x, y,z) [x, '\newline', y,'\newline',z], condition_name1, condition_name2,condition_name3, 'UniformOutput', false);
    set(gca, 'XTickLabel', cat_plot);
    set(gca, 'TickLabelInterpreter', 'tex');

    ylabel('mean RT [ms]');
    ylim([0  ceil(par.max_y/100)*100]);
    RTtitle_name = append('  All Sessions: ', '\newline', ' Reaction Time',...
        '\newline','         dir ',num2str(ndir));
    title(RTtitle_name,'Color','k','Interpreter','tex');
    hold on;
    % Testi da aggiungere per ogni gruppo di errorbar
    texts_e = {'Solo S','Solo K','Joint S','Joint K'};
    for k = 1:length(RT_pldirMean)
        x= k;
        e = errorbar(x, RT_pldirMean(k), RT_pldirSt(k), ...
            'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', colorerror{k}, 'LineWidth', 1);
        e.Annotation.LegendInformation.IconDisplayStyle = 'off';
        if par.idcond==1
            text(x, RT_pldirMean(k) + RT_pldirSt(k) + 5, texts_e{k}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);
        end
    end
    % text(0.95, 0.9,strcat('Chamber:',data_trials(1).Chamber), 'Units', 'normalized', ...
    %     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    %     'FontSize', 12, 'Color','k');
    set(gca, 'FontSize', 14);
    set(gcf, 'WindowState', 'maximized');


    RT_barplot_path = strcat(par.dir_path,'RT_Barplot\Single_dir');

    par.savePlotEpsPdfMat.dir_png = strcat(RT_barplot_path,'\PNGs\');
    par.savePlotEpsPdfMat.dir_pdf = strcat(RT_barplot_path,'\PDFs\');
    par.savePlotEpsPdfMat.dir_mat = strcat(RT_barplot_path,'\MATfiles\');

    name_fig = strcat('RT_barplot_','dir_',num2str(ndir));
    par.savePlotEpsPdfMat.file_name = name_fig;
    savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
end