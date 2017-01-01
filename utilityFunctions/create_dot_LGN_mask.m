function create_dot_LGN_mask(session_dir,func,gfeatNum,zthresh)

% Creates:
%
%   lh.<func>.LGN.nii.gz
%   rh.<func>.LGN.nii.gz
%
%   in the session directory, based on a moving dot functional localizer.
%
%   Usage:
%   create_dot_LGN_mask(session_dir,gfeatNum,zthresh)
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
if ~exist('zthresh','var')
    zthresh = 2.57; % zstat threshold for p < 0.01
end
hemis = {'lh' 'rh'};
%% get bold dirs
d = find_bold(session_dir);
%% Load LGN masks
for hh = 1:length(hemis)
    hemi = hemis{hh};
    LGNmask = load_nifti(fullfile(session_dir,[hemi '.LGN.nii.gz']));
    featDir = fullfile(session_dir,d{gfeatNum},[func '.gfeat']);
    switch hemi
        case 'lh'
            invol = fullfile(featDir,'cope4.feat','stats','zstat1.nii.gz');
        case 'rh'
            invol = fullfile(featDir,'cope3.feat','stats','zstat1.nii.gz');
    end
    % save dotLGN mask
    zstat = load_nifti(invol);
    dotmask = zstat.vol>zthresh;
    goodind = LGNmask.vol > 0 & dotmask;
    dotLGN = LGNmask;
    dotLGN.vol = zeros(size(dotLGN.vol));
    dotLGN.vol(goodind) = 1;
    save_nifti(dotLGN,fullfile(session_dir,[hemi '.' func '.LGN.nii.gz']));
end