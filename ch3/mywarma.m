function [a,gamma]=mywarma(y,n,m,M)
%
% The modified Yule-Walker ARMA method 
% given by equation (3.7.9) with the AR coefficients estimated
% using the overdetermined set of equation (3.7.4), where W is
% the identity matrix.
%
% [a,gamma]=mywarma(y,n,m,M)
% 
%      y     -> the data vector
%      n     -> AR model order
%      m     -> MA model order
%      M     -> the constant which determine the amount of 
%               overdetermination
%      a     <- the AR coefficients vector
%      gamma <- the onesided gamma in equation (3.7.8)

% Copyright 1996 by R. Moses

y=y(:);
N=length(y);             % data length

%-------------------- estimate the AR coefficients --------------------

% compute the standard biased ACS estimate [r(0) r(1) r(2) ...r(n)]
if (N <= m+M)
   disp('Error: m+M must be < the data length.');
   return
elseif (M<n)
   disp('Error: M must be >= n');
   return
end
r=zeros(m+M+1,1);
for i = 0 : m+M,
   r(i+1)=y(1:N-i)'*y(i+1:N)/N;
end

% form the r vector and R matrix in equation (3.7.1)
r1=r(m+2:m+M+1);
if ((m-n+1)>=0)
   R1=toeplitz(r(m+1:m+M),r(m+1:-1:m-n+2).');
else
   R1=toeplitz(r(m+1:m+M),[r(m+1:-1:1);conj(r(2:abs(m-n+1)+1))].');
end

% compute the AR coffecients

a=-R1\r1;
a=[1;a];

%-------------------- estimate the gamma coefficients --------------------
gamma=zeros(m+1,1);
for k = 0 : m,
   for j = 0 : n,
      for p = 0 : n,
         ind=k+p-j;
         if ind >=0
            gamma(k+1)=gamma(k+1)+a(j+1)*a(p+1)'*r(ind+1);
         else
            gamma(k+1)=gamma(k+1)+a(j+1)*a(p+1)'*r(abs(ind)+1)';
         end
      end
   end
end          
