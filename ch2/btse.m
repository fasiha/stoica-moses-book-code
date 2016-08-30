function phi=btse(y,w,L)
%
% The spectral estimator using the Blackman-Tukey method.
% The covariance lags are obtained from the standard biased
% estimate.
%
% phi=btse(y,w,L);
%
%      y -> the data vector
%      w -> the window function (for k=0, ..., M-1)
%      L -> the number of psd samples
%    phi <- spectral estimates at L frequencies w=0, 2*pi/L, ..., 2*pi(L-1)/L

% Copyright 1996 by R. Moses

% make sure the length of the window is less than or equal 
% to the total data length
M=length(w);
N=length(y);
if (M>N)
   error('The length of the window is longer than the data length.');
   return
end

% generate the covariance estimate
r=xcorr(y,'biased');    %vector of biased covariance estimates
r=r(N:N+M-1);           % get r(0) to r(M) only
    
rw = r(:).*w(:);        % the windowed ACS from 0 to M

% generate the spectral estimate

phi=2*real(fft(rw,L))-rw(1);    





