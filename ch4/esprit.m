function w=esprit(y,n,m)
%
% The ESPRIT method for frequency estimation.
%
%  w=esprit(y,n,m);
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
S=U(:,1:n);

phi = S(1:m-1,:)\S(2:m,:);

w=-angle(eig(phi));
return

