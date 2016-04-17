function average_xhemi_vol(session_dir,subject_name,lh_vol_in,rh_vol_in,MNI_out,mh_vol_out,map_type)

% Takes in left and right 'V1corr' files, averages in cvsMNI space, and
% projects back to subject anatomical space
%
%   Usage:
%   average_V1corr(session_dir,subject_name,lh_vol_in,rh_vol_in,MNI_out,mh_vol_out,map_type)
%
%   Written by Andrew S Bock Feb 2016

%% set defaults
if ~exist('map_type','var')
    map_type = 'foo'; % if 'copol' or 'varpol', uses 'circ_mean';
end
%% Average in cvsMNI space
average_in_cvsMNI(session_dir,subject_name,map_type,lh_vol_in,rh_vol_in,MNI_out)
%% Project back to subject space
ref_vol = lh_vol_in;
apply_cvs_inverse(MNI_out,mh_vol_out,ref_vol,subject_name);