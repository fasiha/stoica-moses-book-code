function phi=daniellse(y,J,L)
%
% Spectral estimation using the Daniell method.
%
% phi=daniellse(y,J,L);
%
%      y -> the data vecto
%      J -> 2J+1 = the number of frequency samples to average
%      L -> the number of psd samples
%
%    phi <- spectral estimates at L frequencies w=0, 2*pi/L, ..., 2*pi(L-1)/L

% Copyright 1996 by R. Moses

% compute the peroidogram samples of length L
y1=periodogramse(y,ones(length(y),1),L);

% add elements to the beginning and end of y1 for filter continuity

y1=[y1(L-J+1:L);y1;y1(1:J)];

b=ones(1,2*J+1);      % filter parameters for moving average
a=1;

phi=filter(b,a,y1)/(2*J+1);
phi=phi(2*J+1:length(phi));

