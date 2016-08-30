function phi=capon_sp(Y,L,d)
%
% The Capon method for direction of arrival estimation.
%
% phi=capon_sp(Y,L,d);
%
%    Y   <- the ULA data
%    L   <- the number of samples on [-pi/2,pi/2] to search for the sources
%    d   <- sensor spacing
%    phi -> the spatial spectral estimate at L equally spaced angles
%           in [-90,90] degrees
%

% Copyright 1996 by R. Moses

[m,N]=size(Y);

% compute the inverse of the sample covariance matrix
IR=inv(Y*Y'/N);

phi=zeros(L,1);

for i = 1 : L,
   a=exp(-2*pi*j*d*sin(-pi/2 + pi*(i-1)/L)*[0:m-1].');
   phi(i)=1/real(a'*IR*a);
end
