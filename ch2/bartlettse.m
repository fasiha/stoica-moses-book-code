function phi=bartlettse(y,M,L)
%
% The Bartlett method of spectra estimation.
%
% phi=bartlettse(y,M,L);
%
%      y -> the data vector
%      M -> the length of subsequences of y
%      L -> the number of psd samples
%
%    phi <- spectral estimates at L frequencies w=0, 2*pi/L, ..., 2*pi(L-1)/L

% Copyright 1996 by R. Moses


% check the lenth M
N=length(y);
if (M>N)
   error('M is greater than the data length.');
   return
end

phi=welchse(y,ones(M,1),M,L);   % bartlett is a special case of welch.
