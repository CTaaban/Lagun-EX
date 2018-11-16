function k = wavenumberGuo(h,T,g)
%WAVENUMBERGUO Compute wavenumber according to Guo 2002.
%   K = WAVENUMBERGUO(H,T,g) returns the wave number K, given the total
%   water depth in meters H (which includes wave set-up and tidal
%   influences) the wave period in seconds T, and G (gravitational
%   acceleration of the Earth). T and / or H can be vectors, in which
%   case K is also a vector. This function implements a simple yet
%   accurate approximation formula for the linear dispersion relation.
%   It is based on the work of Guo, as published in Coastal
%   Engineering (2002).

% constants
beta = 2.4908;                                  % Guo's equation 25

sigma = 2 * pi ./ T;                            % Angular frequency
x = h .* sigma ./ sqrt(g .* h);                 % Guo's equation 11, part 1
y = x.^2 .* (1 - exp(-x .^ beta)) .^ (-1/beta); % Guo's equation 21
k = y ./ h;                                     % Guo's equation 11, part 2