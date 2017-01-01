function [opt_HRFparams,fval] = fit_pRF(srctc,TR,X,Y,px,py,sigmaMajor,HRFparams,images,lower_bound,upper_bound)

% Finds the optimal pRF parameters for a given timecourse
%
%   Usage:
%       [opt_params] = fit_pRF(srctc,TR,X,Y,images,init_params,lower_bound,upper_bound)
%
%   Inputs:
%       srctc - source timecourse, TRx1 column vector
%       TR - Repetition time (seconds)
%       X - column vector of X visual angle of stimulus, Nx1 column vector
%       Y - column vector of Y visual angle of stimulus, Nx1 column vector
%       images - NxTR matrix, created using pRF.m
%       init_params - set of initial pRF parameters
%           x0 - center x location of pRF
%           y0 - center y location of pRF
%           sigma - size of pRF
%           peakt - time to peak of HRF
%           undert - time to undershoot of HRF
%       lower_bound - lower bounds on opt_params
%       upper_bound - upper bounds on opt_params
%
%   Output:
%       opt_params - optimal pRF parameters, same convention as init_params
%
%   Written by Andrew S Bock May 2015

%% Set defaults
if ~exist('init_params','var')
    %init_params = [0 0 2 6 16];
    HRFparams = [6 16];
end
if ~exist('lower_bound','var')
    %lower_bound = [-30 -30 0.2 2 4]; % x0,y0,sigma,peakt,undert
    lower_bound = HRFparams - 5; % x0,y0,sigma,peakt,undert
end
if ~exist('upper_bound','var')
    %upper_bound = [30 30 30 10 20]; % x0,y0,sigma,peakt,undert
    upper_bound = HRFparams + 5; % x0,y0,sigma,peakt,undert
end
%% Find the optimal params using GlobalSearch
% problem = createOptimProblem(...
%     'fmincon','objective',...
%     @(params) calc_pRF(srctc,TR,X,Y,params,images), ...
%     'x0',init_params,'lb',lower_bound,'ub',upper_bound,...
%     'nonlcon',@tpcon);
% problem = createOptimProblem(...
%     'fmincon','objective',...
%     @(HRFparams) calc_pRF(srctc,TR,X,Y,px,py,sigmaMajor,HRFparams,images), ...
%     'x0',HRFparams,'lb',lower_bound,'ub',upper_bound);
problem = createOptimProblem(...
    'fmincon','objective',...
    @(HRFparams) calc_pRF(srctc,TR,X,Y,px,py,sigmaMajor,HRFparams,images), ...
    'x0',HRFparams,'lb',lower_bound,'ub',upper_bound,...
    'nonlcon',@tpcon);
gs = GlobalSearch('Display','iter','TolFun',0.001,'TolX',0.001);
[opt_HRFparams,fval] = run(gs,problem);