function [a,sig2]=lsar(y,n)
%
% The Least-Squares AR method (the covariance method)
% given by equation (3.4.14) with N1=n+1 and N2=N.
%
% call [a,sig2]=lsar(y,n);
% 
%      y    -> the data vector
%      n    -> AR model order
%      a    <- the AR coefficient vector estimate
%      sig2 <- the white noise variance estimate

% Copyright 1996 by R. Moses
y=y(:);
N=length(y);             % data length

% compute the standard biased ACS estimate [r(0) r(1) r(2) ...r(n)]
if (N <= n)
   disp('Error: the AR model order is greater than or equal to the data length.');
   return
end
r=zeros(n+1,1);
for i = 0 : n,
   r(i+1)=y(1:N-i)'*y(i+1:N)/N;
end

% form the y vector and Y matrix given in equation (3.4.14)
% with the first and the last n rows removed
y1=[y(n+1:N)];
Y1=toeplitz(y(n:N-1),y(n:-1:1).');

% compute the AR coffecients
a= -Y1\y1;


% compute the noise variance
sig2=norm(Y1*a+y1)^2/(N-n);

a=[1;a];
