function order=sinorder(mvec,mse,N,nu);
% function order=sinorder(mvec,mse,N,nu);
%
% AIC order estimation for sinudoidal model
% 
% inputs:
%  mvec: vector of number of (real or complex) sinusoids
%  mse: vector mean square errors (that is, estimate of sigma^2) for model orders
%       given in mvec.
%  N:   number of real-valued data points 
%  nu: GIC parameter  (usually nu \in [2,6]; default=4

% output:
% order:  the model order that minimizes the AIC criterion.

% randy moses, 09 sep 2003

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
%
% for the first three methods, m=3*n_c+1, where n_c=#sinusoids
% for BIC, m=5*n_c+1, where n_c=#sinusoids

AIC  = N*log(mse) + 2*(3*mvec+1);
AICc = N*log(mse) + (2*(3*mvec+1) * N)./(N-(3*mvec+1)-1);
GIC  = N*log(mse) +  (3*mvec+1) * nu;
BIC  = N*log(mse) +  (5*mvec+1) * log(N);

[a,b]=min(AIC);
order(1)=mvec(b);

[a,b]=min(AICc);
order(2)=mvec(b);

[a,b]=min(GIC);
order(3)=mvec(b);

[a,b]=min(BIC);
order(4)=mvec(b);
return;
