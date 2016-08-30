function phi=correlogramse(y,L)
%
% The correlogram spectral estimator.
%
% phi=correlogramse(y,L)
%
%      y -> the data vector
%      L -> the number of psd samples
%    phi <- spectral estimates at L frequencies w=0, 2*pi/L, ..., 2*pi(L-1)/L

% Copyright 1996 by R. Moses

y=y(:);         % columlize the data matrix
phi=periodogramse(y,ones(size(y)),L);  % the correlogram SE is the same as 
                                       % the unwindowed periodogram SE
