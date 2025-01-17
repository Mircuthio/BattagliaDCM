
function RTbarplotAllvsSingleDir(RT_plotMean,RT_plotStd,RT_Ssingle,RT_KSingle,RTSingle_name,par)

color_name = par.color_name;
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

for ndir=1:par.num_dir
    dirField = strcat('dir', num2str(ndir));
    RT_pldirMean = RT_plotMean(:,ndir);
    RT_pldirSt = RT_plotStd(:,ndir);

    for cd=1:par.num_cond
        RT_singleScond = RT_Ssingle(cd).(strcat('dir',num2str(ndir)));
        RT_singleKcond = RT_KSingle(cd).(strcat('dir',num2str(ndir)));
        for k=1:length(RT_singleKcond)
            figure
            RT_singlePlotS = 1000*RT_singleScond(k);
            RT_singlePlotK = 1000*RT_singleKcond(k);

            RT_finalPlot = [RT_pldirMean;RT_singlePlotS;RT_singlePlotK];
            RT_barplot = bar(RT_finalPlot,'FaceColor', 'flat');
            cond_name1new = {condition_name1{1:end},condition_name1{1:2}};
            if cd == 3
                colorplotnew = {colorplot{1:end},colorplot{3:4}};
                cond_name2new = {condition_name2{1:end},condition_name2{3:4}};
                cond_name3new = {condition_name3{1:end},condition_name3{3:4}};
                texts_eSingle = {'Joint S','Joint K'};
            elseif cd == 1
                colorplotnew = {colorplot{1:end},colorplot{1:2}};
                cond_name2new = {condition_name2{1:end},condition_name2{1},condition_name2{1}};
                cond_name3new = {condition_name3{1:end},condition_name3{1},condition_name3{1}};
                texts_eSingle = {'Solo S','Obs K'};
            else
                colorplotnew = {colorplot{1:end},colorplot{1:2}};
                cond_name2new = {condition_name2{1:end},condition_name2{2},condition_name2{2}};
                cond_name3new = {condition_name3{1:end},condition_name3{2},condition_name3{2}};
                texts_eSingle = {'Obs S','Solo K'};
            end
            for kcol=1:length(RT_finalPlot)
                RT_barplot.CData(kcol,:) = colorplotnew{kcol};
            end
            % infoNaN_plot = infoNaN(cd,:);
            cat_plot = cellfun(@(x, y,z) [x, '\newline', y,'\newline',z], cond_name1new, cond_name2new,cond_name3new, 'UniformOutput', false);
            set(gca, 'XTickLabel', cat_plot);
            set(gca, 'TickLabelInterpreter', 'tex');
            ylabel('RT [ms]');
            ylim([0  ceil(par.max_y/100)*100]);
            RTtitle_name = append('Reaction Time',...
                '\newline','         Dir ',num2str(ndir));
            title(RTtitle_name,'Color',color_name{cd},'Interpreter','tex');
            hold on;
            xline(4.5,'--k')
            % Testi da aggiungere per ogni gruppo di errorbar
            texts_e = {'Solo S','Solo K','Joint S','Joint K'};
            for kk = 1:length(RT_pldirMean)
                x= kk;
                e = errorbar(x, RT_pldirMean(kk), RT_pldirSt(kk), ...
                    'vertical', 'linestyle', 'none', 'CapSize', 5, 'Color', colorerror{kk}, 'LineWidth', 1);
                e.Annotation.LegendInformation.IconDisplayStyle = 'off';
                if par.idcond==1
                    text(x, RT_pldirMean(kk) + RT_pldirSt(kk) + 5, texts_e{kk}, ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);
                end
            end
            for kke = 5:6
                x = kke;
                if par.idcond==1
                    text(x, RT_finalPlot(kke) + 5, texts_eSingle{kke-4}, ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);
                end
            end
            % text(0.95, 0.9,strcat('Chamber:',data_trials(1).Chamber), 'Units', 'normalized', ...
            %     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            %     'FontSize', 12, 'Color','k');
            set(gca, 'FontSize', 14);
            set(gcf, 'WindowState', 'maximized');
            textString = sprintf("Session: %s \nChamber: %s\nCondition: %d",par.session_name,par.chamber, cd);

            textString = sprintf("%s\nD%d-T%d", textString, RTSingle_name(cd).(dirField).D, RTSingle_name(cd).(dirField).T(k));
            
            textString2 = sprintf("Session: All \nChamber: All\nCondition: All \nD%d-T All",ndir);
            text(0.15, 0.95, textString2,  'Units', 'normalized',...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                'FontSize', 12, 'Color', 'k');

            text(0.95, 0.95, textString,  'Units', 'normalized',...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                'FontSize', 12, 'Color', 'k');
            
            RT_barplot_path = strcat(par.dir_path,'RT_Barplot\ConfrontoTrials');

            par.savePlotEpsPdfMat.dir_png = strcat(RT_barplot_path,'\PNGs\');
            par.savePlotEpsPdfMat.dir_pdf = strcat(RT_barplot_path,'\PDFs\');
            par.savePlotEpsPdfMat.dir_mat = strcat(RT_barplot_path,'\MATfiles\');

            name_fig = strcat('RT_barplot_','dir_',num2str(ndir),'_Trial_',num2str(RTSingle_name(cd).(dirField).T(k)));
            par.savePlotEpsPdfMat.file_name = name_fig;
            savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
            close all
        end
    end
end