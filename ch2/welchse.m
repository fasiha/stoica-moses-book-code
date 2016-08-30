function phi=welchse(y,v,K,L)
%
% The Welch method of spectral estimation.
%
% phi=welchse(y,v,K,L);
%
%      y -> the data vector 
%      v -> the window vector
%      K -> (j-1)K+1 is the starting point for the jth subsequence
%      L -> the number of psd samples
%
%    phi <- spectral estimates at L frequencies w=0, 2*pi/L, ..., 2*pi(L-1)/L

% Copyright 1996 by R. Moses

N=length(y);     % total data length
M=length(v);     % length of each observations to split from y

% check the lengths
if (M > N)
   error('The window length is larger than the data length.');
   return
elseif (K > N)
   error('The value of K exceeds the data length.');
   return
end

S=fix((N-M+K)/K); % number of subsamples
P=mean(v.^2);     % the power of the window vector v

y=y(:);          % make y a column vector

% compute the weighted periodgram for each subsample of observations
% and sum the results

phi=zeros(L,1);
for i = 1 : S,
   phi=phi + periodogramse(y((i-1)*K+1:(i-1)*K+M),v,L);
end

phi=phi/S/P;
