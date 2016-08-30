function freqplot(w,a,wh,ah)
% Plots line spectra for real-valued signals given the frequencies and complex amplitudes.
% 
% Inputs:
%   w, a: real-valued positive frequencies and corresponding amplitudes of the true signal.
%   wh, ah: real-valued freqency estimates (positive and negative pairs) and corresponding
%          complex-valued amplitude estimates.
%   w and wh are in radians, and take on values in (0,pi) for w and
%   [-pi,pi] for wh
% Output:
%  a line plot on the current figure in which the true frequencies and
%  amplitudes are shown as dotted vertical lines, and estimates are shown
%  as solid vertical lines.

% first the wh and ah 
wh=wh(:);
ah=ah(:);

ind=find(wh>=0 & wh<=pi);
fh=wh(ind)/pi/2;
ah=2*abs(ah(ind));

fh=fh(:)';
ah=ah(:)';
x=[fh;fh]; y=[0*ah;ah];
line(x,y,'Color','k')

% now overlay the dotted lines for w and a 

f=w(:)'/pi/2;
a=a(:)';
x=[f;f]; y=[0*a;a];
line(x,y,'Color','r','Linestyle',':')

axis([0,.5,0,1.1*max(max(ah),max(a))]);



