function [a,sig2]=yulewalker(y,n)
% 
% The Yule-Walker method for AR spectral estimation, given
% by equation (3.4.2).
%
% [a,sig2]=yulewalker(y,n);
% 
%      y    -> the data vector
%      n    -> AR model order
%      a    <- the AR coefficient vector estimate
%      sig2 <- the white noise variance estimate

% Copyright 1996 by R. Moses

y=y(:);  
N=length(y);             % data length

if (N < n)
   disp('Error: the AR model order is greater than the data length.');
   return
end

% compute the standard biased ACS estimate [r(0) r(1) r(2) ...r(n)]
%
r=zeros(n+1,1);
for i = 0 : n,
   r(i+1)=y(1:N-i)'*y(i+1:N)/N;
end

% form the Toeplitz covariance matrix
Rn=toeplitz(conj(r(1:n)));

% compute the AR coffecients
a=-Rn\r(2:n+1);

% compute the noise variance
sig2=real(r(1)+a.'*conj(r(2:n+1)));

a=[1;a];
