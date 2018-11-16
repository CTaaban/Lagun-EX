function [H13] = significant_wave_height(Hind)
% [H13] = significant_wave_height(Hind)
% Compute the signicant wave height 
% input: Hind is vector containing the individual wave heights of a given record (m) 
% output: H13 significant wave height (m)

n_waves = length(Hind); % number of waves in the signal

% 1. sort the wave heights in descending order
H_sort = sort(Hind,'descend'); 

% 2. define a new vector H_new containing the third highest waves of the record
% (=first n13 elements of the sorted vector)
n13 = round(n_waves/3); % the round function insures that n13 is an integer 
H_new = H_sort(1:n13);

% 3. calculate H13 (output of the function)
H13 = sum(H_new)/n13;