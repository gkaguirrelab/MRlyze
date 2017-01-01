function [neg_co,PredResp] = calc_pRF_pred(srctc,TR,X,Y,px,py,psig,HRFparams,images)

%% Pull out params
peakt = HRFparams(1);
undert = peakt+10; %HRFparams(2);
%% Make HRF
hrf =  doubleGammaHrf(TR,[peakt undert]);
%% Translate grid so that center is at RF center
nX = X - px;   % positive x0 moves center right
nY = Y - py;   % positive y0 moves center up
%% make gaussian on current grid
%cG =   sigs(:,:,3).*(exp (-(nY.^2 + nX.^2) ./ (2*sigs(:,:,1).^2)));
%sG =   sigs(:,:,4).*(exp (-(nY.^2 + nX.^2) ./ (2*sigs(:,:,2).^2)));
%rf = cG - sG;
rf = exp (-(nY.^2 + nX.^2) ./ (2*psig(1).^2));
%% now convolve with Hrf
allstimimages = filter(hrf,1, images');
%% convolve with stimulus
PredResp = allstimimages*rf;
% Set timecourses with very little variation (var<0.1) to flat
%PredResp = set_to_flat(PredResp);
%% Compute negative variance explained
%   We want the maximize the correlation, but as we will use
%   fminsearch later (which finds the minimum) we convert to a negative value.
if ~isempty(srctc)
    co = corr2(srctc,PredResp);
    neg_co = -co;
else
    neg_co = nan;
end
