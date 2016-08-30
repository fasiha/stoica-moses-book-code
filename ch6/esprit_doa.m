function doa=esprit_doa(Y,n,d)
%
% The ESPRIT method for DOA estimates 
%
% call doa=esprit_doa(Y,n,d)
%
%      Y    <- the ULA data
%      n    <- the number of sources
%      d    <- sensor spacing in wavelengths
%      doa  -> the vector of DOA estimates

% Copyright 1996 by R. Moses

[m,N]=size(Y);

% compute the sample covariance matrix
R=Y*Y'/N;

% do the eigendecomposition; use svd because it sorts eigenvalues
[U,D,V]=svd(R);

S=U(:,1:n);

phi = S(1:m-1,:)\S(2:m,:);
w=-angle(eig(phi));
doa=asin(w/d/pi/2)*180/pi;
