
function [c,ceq] = sigcon(x)

% Used in calc_CF. Contrains the DoG sigmas, such that the negative
% Gaussian sigma is always larger than the positive Gaussian sigma
%
%   Written by Andrew S Bock Apr 2015

% Sigma 1 > Sigma 2 (center-surround)
CS = x(1) - x(2);

% All nonlinear contraints
c = [CS];

% No equality contraints
ceq = [];

% % Sigma 1 > Sigma 2 (center-surround)
% CS = x(1) - x(2);
% 
% % Contrain the vertex index (sig(4)) to be an integer
% integer = x(4);
% 
% % All nonlinear contraints
% c = [CS;integer];
% 
% % No equality contraints
% ceq = [];