function phi = argamse(gamma,a,L);  
% generates L samples of an ARMA spectral density function phi
% from the ARMA coefficients.
%
% phi = argamse(gamma,a,L);  
%     gamma -> the spectral density numerator coefficient vector 
%              [gamma(0),..., gamma(m)]^T
%     a     -> the AR coefficient vector (including the leading '1')
%     phi   <- the spectral density at frequencies 0, 2pi/L, ... 2pi*(L-1)/L

% Copyright 1996 by R. Moses
a=a(:);  gamma=gamma(:);
m=length(gamma)-1;
H=freqz(1,a,L,'whole');
num = real(fft([gamma;zeros(L-2*m-1,1);conj(gamma(m+1:-1:2))]));
phi = num.*(abs(H).^2);
