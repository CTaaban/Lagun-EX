function [E, f, confLow, confUpper] = wave_spectrum(data,nfft,Fs)
% [E, f, confLow, confUpper] = wave_spectrum(data,nfft,Fs)
% Compute variance spectral density spectrum of the time-series and its 
% 90% confidence intervals. 
% The time-series is first divided into blocks of length nfft before being 
% Fourier-transformed.
%
% INPUT
%   data    timeseries 
%   nfft    block length
%   Fs      sampling frequency (Hz)
%
% OUTPUT
%   E  variance spectral density. If data is in meters, E is in m^2/Hz
%   f  frequency axis (Hz)
%   confLow and confUpper  Lower and upper 90% confidence interval; 
%                          (Multiplication factors for E)


% 1. PRELIMINARY CALCULATIONS
% ---------------------------
data = data(:);
n = size(data,1);            % length of the time-series
nfft = nfft - rem(nfft,2);   % the length of the window should be an even number

data = detrend(data);        % detrend the time-series
nBlocks = floor(n/nfft);     % number of blocks (use of floor to make it an integer)

data_new = data(1:nBlocks*nfft); % (we work only with the blocks which are complete)


% we organize the initial time-series into blocks of length nfft 
dataBlock = reshape(data_new,nfft,nBlocks); % each column of dataBlock is one block


% 2. CALCULATION VARIANCE DENSITY SPECTRUM
% ----------------------------------------

% definition frequency axis
df = Fs/nfft; % frequency resolution of the spectrum df = 1/[Duration of one block]
f = (0:df:Fs/2)'; % frequency axis (Fs/2 = Fnyquist = max frequency)
fId = (1:length(f))'; 

% Calculate the variance for each block and for each frequency
fft_data = fft(dataBlock,nfft,1);     % Fourier transform of the data 
fft_data = fft_data(fId,:);           % Only one side needed
A = 2/nfft*real(fft_data);            % A(i,b) and B(i,b) contain the Fourier coefficients Ai and Bi for block b 
B = 2/nfft*imag(fft_data);            %  -- see LH's book, page 325 for definition of Ai and Bi
                                      % /!\ assumes that mean(eta)=0 

E = (A.^2 + B.^2)/2;                  % E(i,b) = ai^2/2 = variance at frequency fi for block b. 

% We finally average the variance over the blocks, and divide by df to get the variance DENSITY spectrum
E = mean(E,2)/df;              

% 3. CONFIDENCE INTERVALS
% -----------------------
edf = round(nBlocks*2); % Degrees of freedom 
alpha = 0.1; % calculation of the 90% confidence interval
confLow = edf/chi2inv(1-alpha/2,edf); % see explanations on confidence intervals given in lecture 3 
confUpper  = edf/chi2inv(alpha/2,edf);
