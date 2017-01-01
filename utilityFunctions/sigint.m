function [c,ceq] = sigint(x)

% Used in calc_CF. Contrains the DoG sigmas, such that the positive
% Gaussian sigma is always larger than the negative Gaussian sigma
%
%   Written by Andrew S Bock Apr 2015

c = isint(x(4)); 
ceq = [];