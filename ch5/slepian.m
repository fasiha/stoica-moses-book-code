function h=slepian(N,K,J)
%
% Generates the first J Slepian sequences of length N and "bandwidth 
% factor" K as defined in Section 5.3.
%
% h=slepian(N,K,J);
%
% N <- the length of each Slepian sequence
% K <- the ratio of the baseband filter bandwidth to 1/N
% J <- the number of the Slepian sequences to return
% h -> an NxJ matrix whose ith column is the ith Slepian sequence
%

% Copyright 1996 by R. Moses

% generate the Gamma matrix 
Gamma = toeplitz(K/N*sinc(K/N*[0:N-1]));

% eigenvalue decomposition; use svd because it orders eigenvalues
[U,D,V]=svd(Gamma);

h=U(:,1:J);

% this next line is optional, but ensures that the first 
% Slepian sequence has a positive center element (as is usually
% the case for a lowpass filter).
if sum(h(:,1))<0, h(:,1)=-h(:,1); end
