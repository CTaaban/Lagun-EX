function Hrms = rms_wave_height(Hind)
% Hrms = rms_wave_height(Hind)
% input    Hind array containing the individual wave heights in m 
% output   Hrms root mean square wave height in m.

Hrms = sqrt(mean(Hind.^2)); % calculation of Hrms (NB: sqrt is a matlab built-in function that calculate the root square)
