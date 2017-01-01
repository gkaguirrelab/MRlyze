function [hrf] = doubleGammaHrf(TR,tp,beta,rt,len,scale)

% Computes a double gamma Hemodynamic Response Function (HRF).
%
%   Usage: [hrf] = doubleGammaHrf(TR,tp,fwhm,alp,len)
%
%   Defaults:
%       TR = 2; Repetition time
%       tp = [5 15]; time to peak/undershoot
%       beta = [1 1]; scale of peak/undershoot
%       rt = 1/10; ratio of response to undershoot
%       len = 33; length of HRF (in seconds)%
%
%   Equations modeled after:
%       DeSimone, Viviano, Schneider (2015) J Neuro.
%
%   Written by Andrew S Bock May 2015

%% Set defaults
if ~exist('TR','var')
    % Repetition time
    TR = 1;
end
if ~exist('tp','var')
    % time to peak response
    tp(1) = 5;
    % time to minimum undershoot
    tp(2) = 15;
end
if ~exist('beta','var');
    % peak
    beta(1) = 1;
    % undershoot
    beta(2) = 1;
end
if ~exist('rt','var')
    % ratio of response to undershoot
    rt = 1/10;
end
if ~exist('len','var')
    % length of HRF
    len = 33;
end
if ~exist('scale','var')
    scale = 1;
end
%% Create HRF
alpha_1 = tp(1)/TR;
alpha_2 = tp(2)/TR;
beta_1 = beta(1);
beta_2 = beta(2);
t = TR:TR:len;
hrf = scale*...
    (...
    ((t.^(alpha_1) * beta_1 ^ alpha_1 .* exp(-beta_1 * t)) / ...
    gamma(alpha_1)) ...
    - ...
    rt * ((t.^(alpha_2) * beta_2 ^ alpha_2 .* exp(-beta_2*t))/ ...
    gamma(alpha_2)) ...
    );
hrf = hrf / trapz(hrf);