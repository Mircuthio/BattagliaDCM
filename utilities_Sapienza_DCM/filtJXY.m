%% Filtraggio passa-basso per dati JXY

function data_out = filtJXY(data,par)

fs = 1000;
Nc = par.Nc;
Fc = par.Fc;

LP = fir1(Nc, Fc/(fs/2));
data_out = filtfilt(LP, 1,data);
