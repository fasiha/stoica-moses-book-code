function order=orderest(mo,sig2,N,nu);
% function orderest(mo,sig2,N,nu);
%
% Order estimation for a generic ARMA model
% 
% inputs:
%  mo: vector of model orders
%  sig2: vector mean square errors (that is, estimate of sigma^2) for model orders
%       given in mvec.
%  N:   number of data points 
%  nu: GIC parameter  (usually nu \in [2,6]; default=4

% output:
% order:  the model order that minimizes the order estimation criterion.

% randy moses, 05 may 2005

if nargin<4, 
    nu=4;
end

% the general order estimation rule is -2 ln p_m + eta(m,N)*m
%
% where -2 ln p_m = N*ln(sigma^2_m) and where
%
%  for AIC:      eta(m,N) = 2 * m
%  for AIC_c:    eta(m,N) = 2 * (N)/(N-m-1) * m
%  for GIC:      eta(m,N) = nu * m 
%  for BIC(MDL): eta(m,N) = ln(N)* m

AIC  = N*log(sig2) + 2*mo;
AICc = N*log(sig2) + (2*mo * N)./(N-mo-1);
GIC  = N*log(sig2) +  mo * nu;
BIC  = N*log(sig2) +  mo * log(N);

[a,b]=min(AIC);
order(1)=mo(b);

[a,b]=min(AICc);
order(2)=mo(b);

[a,b]=min(GIC);
order(3)=mo(b);

[a,b]=min(BIC);
order(4)=mo(b);

return;
