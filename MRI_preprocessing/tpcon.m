
function [c,ceq] = tpcon(x)

% Used in calc_pRF. Contrains the time to peak to be shorter than the time
% to undershoot by at least 1s.
%
%   Written by Andrew S Bock May 2015

% Time to peak < Time to undershoot (center-surround)
CS = x(1) - x(2) + 1;

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