function doa=root_music_doa(Y,n,d)
%
% The MUSIC method for direction of arrival estimation
%
% call doa=root_music_doa(Y,n,d)
%
%      Y    <- the ULA data
%      n    <- the number of sources
%      d    <- sensor spacing in wavelengths
%      doa  -> the vector of DOA estimates

% Copyright 1996 by R. Moses; last update 7/8/2005

[m,N]=size(Y);

% compute the sample covariance matrix
R=Y*Y'/N;

% do the eigendecomposition; use svd because it sorts eigenvalues
[U,D,V]=svd(R);

G=U(:,n+1:m);

C=G*G';   % mxm
% find the coefficients of the polynomial in (4.5.16)
for kk=1:2*m-1,
    a(kk,1)=sum(diag(C,kk-m));
end
ra=roots(a);

% find the n roots of the a polynomial that are nearest and inside the unit circle,

[dum,ind]=sort(abs(ra));
rb=ra(ind(1:m-1));

% pick the n roots that are closest to the unit circle
[dumm,I]=sort(abs(abs(rb)-1));
w=angle(rb(I(1:n)));


% compute the doas
doa=asin(w/d/pi/2)*180/pi;


return
