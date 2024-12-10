%% function for boxplot execution

function Out_par = boxplotparameters(data,campion,sigma_parameter)

%% parameters evaluation
data = data(~isnan(data));

RT_mean = mean(data);

if campion == 1
    % Standard deviation c
    RT_dev = std(data)/sqrt(length(data));
else 
    RT_dev = std(data);
end
% median
RT_median = median(data);

% max and min
RT_max = max(data);


RT_min = min(data);

% quartili
percentiles = [0.25, 0.75];

RT_quantiles = quantile(data, percentiles);



%% wiskers group of median
lower_whisker_RT_quantile = min(data(data >= RT_quantiles(1) - 1.5*(RT_quantiles(2)-RT_quantiles(1))));
upper_whisker_RT_quantile = max(data(data <= RT_quantiles(2) + 1.5*(RT_quantiles(2)-RT_quantiles(1))));


%% wiskers group of mean with min and max
lower_whisker_RT_mean = (RT_mean - 3*sigma_parameter*RT_dev);
upper_whisker_RT_mean = (RT_mean + 3*sigma_parameter*RT_dev);


% %% wiskers group of mean with min and max
% lower_whisker_RT_mean = min(data(data >= (RT_mean - sigma_parameter*RT_dev) - 2*(2*sigma_parameter*RT_dev)));
% upper_whisker_RT_mean = max(data(data <= (RT_mean + sigma_parameter*RT_dev) + 2*(2*sigma_parameter*RT_dev)));


Out_par.mean = RT_mean;
Out_par.dev = RT_dev;
Out_par.max = RT_max;
Out_par.min = RT_min;
Out_par.median = RT_median;
Out_par.quantile = RT_quantiles;
Out_par.lower_wisker_quantile = lower_whisker_RT_quantile;
Out_par.upper_wisker_quantile = upper_whisker_RT_quantile;
Out_par.lower_wisker_mean= lower_whisker_RT_mean;
Out_par.upper_wisker_mean= upper_whisker_RT_mean;
