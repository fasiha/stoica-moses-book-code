function w=hoyw(y,n,L,M)
%
% The Higher-Order Yule-Walker method.
%
%  w=hoyw(y,n,L,M);
%
%   y  ->  the data vector
%   n  ->  the AR model order
%   L  ->  the number of columns in the matrix in equation (4.4.8)
%   M  ->  the number of rows in the matrix in equation (4.4.8)
%   w  <-  the frequency 

% Copyright 1996 by R. Moses

y=y(:); 
N=length(y);                       % data length
%y1=y(:)-mean(y);

% generate the onesided covariance estimate [r(1),..,r(L+M)]^T
r=zeros(L+M,1);
for i = 1 : L+M,
   r(i)=y(1:N-i)'*y(i+1:N)/N;
end

% form the covariance matrix in (4.4.8)
Omega= toeplitz(r(L:L+M-1),r(L:-1:1).') ;

[U,S,V]=svd(Omega);

% find the matrix for the approximated Omega with rank n
U1=U(:,1:n);
S1=S(1:n,1:n);
V1=V(:,1:n);

% compute the estimate of the b polynomial
b = -V1 * inv(S1) * U1' * r(L+1:L+M);

% find the n roots of the a polynomial that are nearest the unit circle,
% and compute the their angular positions 
rb=roots([1; b]);
[dumm,I]=sort(abs(abs(rb)-1));
w=angle(rb(I(1:n)));

