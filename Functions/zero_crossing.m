function [H_ind , T_ind]=zero_crossing(data,Fs)
% [H_ind , T_ind]=zero_crossing(data,Fs)
% Zero crossing analysis of wave data
% inputs   data: detrended time-series of water elevation in m. 
%          Fs: sampling frequency of data in Hz
% outputs  H_ind: individual wave heights in m (array)
%          T_ind: individual wave periods in s (array)

% 1. Preliminary calculations
% ---------------------------

% time vector for the water elevation [s]
time = (0:length(data)-1)/Fs; 
% Before performing the analysis, we remove the zero values
d0 =data(data~=0);
t0 = time(data~=0);

% 2. Zero-crossing analysis 
% -------------------------

% We identify the indices at which surface elevation changes sign
crossing=find(d0(1:end-1).*d0(2:end)<0);

% If the elevation is negative at t=0, the first crossing
% is a zero-upward-crossing -> it is rejected
if d0(1)<0
    crossing(1)=[];
end
crossing=crossing(1:2:end);  % these are the zero-down-crossing

% 3. Calculation individual wave characteristics
% ----------------------------------------------

if length(crossing)>=2 % Calculate wave period and height if at least one wave has been identified (=2 crossings)
    for i=1:length(crossing)-1
        elevation_crest(i)= max(d0(crossing(i):crossing(i+1)));    % crest elevation
        elevation_trough(i)= min(d0(crossing(i):crossing(i+1)));   % trough elevation
    end  
    T_ind=diff(t0(crossing)); % period = time difference between two successive down-zero-crossing
    H_ind= elevation_crest - elevation_trough; % individual wave heigh
else % if no waves, returns empty vectors
    T_ind =[];
    H_ind=[];
end


