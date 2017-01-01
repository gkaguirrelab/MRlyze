function [outVol] = compare_templates(tName1,tName2,type,tName3,tName4)

% Compares two retinotopic templates
%
%   Written by Andrew S Bock Oct 2015

%% set defaults
if ~exist('type','var')
    type = 'ecc';
    disp('Assuming template is ''ecc''');
end
%% Load in templates
t1 = load_nifti(tName1);
t2 = load_nifti(tName2);
if strcmp(type,'dist')
    t3 = load_nifti(tName3);
    t4 = load_nifti(tName4);
end
%% Compute differences
if strcmp(type,'pol')
    outVol = abs(circ_dist(t1.vol,t2.vol));
elseif strcmp(type,'dist')
    % assume ecc are 1 and 2, pol are 3 and 4
    [x1,y1] = pol2cart(t3.vol,t1.vol);
    [x2,y2] = pol2cart(t4.vol,t2.vol);
    outVol = sqrt( (x1 - x2).^2 + (y1 - y2).^2 );
else
    outVol = abs(t1.vol - t2.vol);
end
disp('done.');
