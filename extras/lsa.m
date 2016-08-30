function [a,rel_MSE]=lsa(y,w)
%
% lsa.m computes the complex amplitudes of 
% sinusoidal components given the frequencies by a least-squares
% solution to the data.
%
% Model:
%   y(n)=sum_{k=1}^p  a_k e^{i w(k) n} + e(n),   n= 0,...,N-1 
%
% Inputs:
%   y   -> Nx1 data vector
%   w   -> sinusoidal frequencies
% Outputs: 
%   a   <- the complex amplitude 
%   rel_MSE <- the relative error between the data and its reconstruction
%   using frequencies and amplitudes

y=y(:);
w=w(:);

A=exp(sqrt(-1)* (0:length(y)-1)'*w');
a=A\y;

yr=A*a;
rel_MSE=1-yr'*yr/(y'*y);

return;