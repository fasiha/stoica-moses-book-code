function phi = armase(b,a,sig2,L);
% generates L samples of an ARMA spectral density function phi
% from the ARMA coefficients.
%
% phi = armase(b,a,sig2,L);
%      b    -> the MA coefficient vector (including the leading '1')
%      a    -> the AR coefficient vector (including the leading '1')
%      sig2 -> the white noise variance 
%     L     -> the number of PSD samples desired
%      phi  <- the spectral density at frequencies 0, 2pi/L, ... 2pi*(L-1)/L

% Copyright 1996 by R. Moses

H=freqz(b,a,L,'whole');
phi = abs(H).^2*sig2;
