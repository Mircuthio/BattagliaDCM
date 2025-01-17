

function maximum = maxdistance(data1_x,data1_y,data2_x,data2_y,tg_pos,tg_pos_names)

y_max = NaN(1,length(tg_pos));
for ntarg = 1:length(tg_pos)
    x_target = tg_pos(ntarg, 1);
            y_target = tg_pos(ntarg, 2);
            if strcmp(tg_pos_names(ntarg),'C')
                % Calcolo della distanza rispetto al centro (0,0) per MS e MK
                data1_distance = sqrt(data1_x.^2 + data1_y.^2);
                data2_distance = sqrt(data2_x.^2 + data2_y.^2);
                y_max(ntarg) = max(max(data1_distance),max(data2_distance));
            else
                data1_distance = sqrt((data1_x - x_target).^2 + (data1_y - y_target).^2);
                data2_distance = sqrt((data2_x - x_target).^2 + (data2_y - y_target).^2);
                y_max(ntarg) = max(max(data1_distance),max(data2_distance));
            end
end
maximum = max(y_max);