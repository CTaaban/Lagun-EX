function [mn] = spectral_moment(f,E,f1,f2,n)
% [mn] = spectral_moment(f,E,fmin,fmax,n)
% calculate the n^th-order spectral moment for a given frequency band [fmin fmax]
% inputs
%  E: variance density spectrum
%  f: frequency axis
%  f1 and f2: minimum and maximum frequency considered in the moment
%  calculation
%  n: order of the moment (if n=0, mn=m0=variance)
% output
%  mn: spectral moment
 
 
% calculation of the spectral moment
if n>=0
    % indices of the frequencies larger than fmin and smaller that fmax
    ind_f = (f>=f1)&(f<=f2);
 else % when n<0, f cannot be equal to zero (as f^(-N)=(1/f)^(N) = infinity if f=0!)
    ind_f = (f>=f1)&(f<=f2)&(f~=0);
 end
 
 mn = trapz(f(ind_f),E(ind_f).*f(ind_f).^n);