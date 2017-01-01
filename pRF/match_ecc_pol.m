function [ind,ecc,pol] = match_ecc_pol(ObsEcc,ObsPol,TrgEcc,TrgPol)

% Finds the closest matching eccentricty and polar angle vertex (Target)
%   within an ROI, given a single eccentricity and polar angle (Observed)
%
%   Written by Andrew S Bock Jul 2015

%% Find the closest matching ecc and pol pair
[ObsX,ObsY] = pol2cart(ObsPol,ObsEcc);
[TrgX,TrgY] = pol2cart(TrgPol,TrgEcc);
diffs = sqrt((ObsX - TrgX).^2 + (ObsY - TrgY).^2);
[~,ind] = min(diffs);
ecc = TrgEcc(ind);
pol = TrgPol(ind);
%%
% V1eccdiff = abs(TrgEcc - ObsEcc);
% V1poldiff = abs(TrgPol - ObsPol);
% [~,eccI] = sort(V1eccdiff);
% [~,polI] = sort(V1poldiff);
% ranklist = 1:length(eccI);
% eccrank(eccI) = ranklist;
% polrank(polI) = ranklist;
% sumrank = eccrank + polrank;
% [~,ind] = min(sumrank);
% ecc = TrgEcc(ind);
% pol = TrgPol(ind);