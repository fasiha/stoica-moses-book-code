function w=music(y,n,m)
%
% The Root MUSIC method for frequency estimation.
%
%  w=music(y,n,m);
%
%      y  ->  the data vector
%      n  ->  the model order
%      m  ->  the order of the covariance matrix in (4.5.14)
%      w  <-  the frequency estimates
%

% Copyright 1996 by R. Moses

y=y(:);
N=length(y);                       % data length

% compute the sample covariance matrix

R=zeros(m,m);
for i = m : N,
   R=R+y(i:-1:i-m+1)*y(i:-1:i-m+1)'/N;
end

% to use the forward-backward approach, uncomment the next line
% R=(R+fliplr(eye(m))*R.'*fliplr(eye(m)))/2;

% get the eigendecomposition of R; use svd because it sorts eigenvalues
[U,D,V]=svd(R);
G=U(:,n+1:m);

C=G*G';   % mxm
% find the coefficients of the polynomial in (4.5.16)
for kk=1:2*m-1,
    a(kk,1)=sum(diag(C,kk-m));
end
ra=roots(a);

% find the m-1 roots of the a polynomial that are nearest and inside the unit circle,
% occasionally, due to numerical inaccuracies, you will get two roots on the
% unit circle and at slightly different angles, instead of two roots at the same
% angle and at slightly different amplitudes.  For this reason we first sort the root
% magnitudes and take the first m-1, then take the n closest of these

[dum,ind]=sort(abs(ra));
rb=ra(ind(1:m-1));

% pick the n roots that are closest to the unit circle
[dumm,I]=sort(abs(abs(rb)-1));
w=angle(rb(I(1:n)));
return

