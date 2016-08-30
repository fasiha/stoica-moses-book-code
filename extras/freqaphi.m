function [a,phi]=freqaphi(y,w)
%
% Finds the least-squares amplitude and phases for sinusoidal data
%  once the frequencies have been estimated.  Uses equation (4.3.8)
%
%  [a,phi]=freqaphi(y,w)
%
%      y   ->  the data vector
%      w   ->  the frequency estimates
%      a   <-  the amplitude estimates
%      phi <-  the phase estimates
%

% Copyright 1996 by R. Moses

y=y(:);
N=length(y);                       % data length

w=w(:); n=length(w);

B= (ones(N,1)*(exp(j*w.'))) .^( [0:N-1]'*ones(1,n) );

beta=B\y;

a  = abs(beta);
phi = angle(beta);

