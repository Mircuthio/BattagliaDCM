function Trials = TR_filterSignals(Trials,CFc,CFe,Nc,Ne)
% Trials = TR_filterSignals(Trials,CFc,CFe,Nc,Ne)
%
% low-pass filter all behavioral signals in all trials
% cut-frequencies for cursor (CFc) and eye (CFe) signals can be set independently
% if only one cut-frequency is entered it will be used for both signals (eye and cursor)
% if no cut-frq is defined, 30Hz and 60Hz will be used for the joystick and
% eye signal respectively (both with a filter order of 100)
%
% NOTE on the choice of the filter order (Nc, Ne):
% ! pay attention when using 'very low' cut frequencies (then N must be 'high')
% To compute the best filter order based on a specific filter design, use one of the following tools:
% - online FIR filter design: http://t-filter.appspot.com/fir/index.html
% - matlab FIR filter design tool: http://www.mathworks.co.uk/help/signal/ug/fir-filter-design.html
%
% sample usage:
% CFc = 25; Nc = 100;
% CFe = 100; Ne = 100;
% TrialsLP = TR_filterSignal(Trials,CFc,CFe,Nc,Ne);
%
% v01  20140124
% v02  20140205  changed order of input parameters
% v03  20140603  TR
% v04  20140604  if one trial's data is too short -> empty JXYEXY

if nargin==1,
    CFc = 30; CFe = 60; Nc = 100; Ne = Nc;
elseif nargin==2,
    CFe = CFc; Nc = 100; Ne = Nc;
elseif nargin==3,
    Nc = 100; Ne = Nc;
elseif nargin==4,
    Ne = Nc;
end


Fs = 1000;  % sampling frequency

% initialize FIR filter
if exist('filtfilt','file')>0,
    LPc = fir1(Nc,CFc/(Fs/2));
    LPe = fir1(Ne,CFe/(Fs/2));
else
    error('Filter function not found - signal processing toolbox not installed?');
end

nTrials = length(Trials);
for tr=1:nTrials,
    
    % if not enough data is available, remove behavioral data (useless)
    L = size(Trials(tr).JXYEXY,2);
    if L<=Nc*3 || L<=Ne*3,
        Trials(tr).JXYEXY = [];
        disp(['--WARNING-- Not enough data points to apply FIR filter in trial ',num2str(tr),'. L = ',num2str(L),', Nc = ',num2str(Nc),', Ne = ',num2str(Ne)]);
        continue
    end
    
    
    % cursor signals
    for s=[1 2 5 6],
        sig = double(Trials(tr).JXYEXY(s,:));
        if length(sig)>Nc,
            % tr
            % s
            if ~any(isnan(sig))
            Trials(tr).JXYEXY(s,:) = filtfilt(LPc,1,sig);
            end
        else
        end
    end
    
    % eye signals
    for s=[3 4 7 8],
        sig = double(Trials(tr).JXYEXY(s,:));
            Trials(tr).JXYEXY(s,:) = filtfilt(LPe,1,sig);
    end
    
end


end