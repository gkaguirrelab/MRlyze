function prepare_pRF_Mathematica(session_dir,subject_name,outName)

% Converts .nii.gz pRF cortical maps that result from 'average_pRF.m' to
% .mgh format, as well as changes the naming convention to that used by
% Mathematica pRF template fitting notebook
%
%   Usage:
%   prepare_pRF_Mathematica(session_dir,subject_name,outName)
%
%   Example:
%   session_dir = '/data/jet/abock/data/Template_Retinotopy/AEK/10012014';
%   subject_name = 'AEK_09242014_MPRAGE_ACPC_7T';
%   outName = 'AEK';
%   prepare_pRF_Mathematica(session_dir,subject_name,outName);
%
%   Written by Andrew S Bock Jun 2015

%% Set defaults
if ~exist('sym','var')
    sym = 1;
end
hemis = {'lh' 'rh' 'mh'};
inmaps = {'co' 'copol' 'coecc'};
outmaps = {'co' 'pol' 'ecc'};
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
%% Convert maps
disp('Converting maps...');
for m = 1:length(inmaps);
    for hh = 1:length(hemis)
        for sy = 1:2
            if sy == 1 % prepare the data on the fsaverage_sym surface
                infile = fullfile(session_dir,[hemis{hh} '.cortex.avg.' inmaps{m} '.prfs.sym.nii.gz']);
                outfile = fullfile(SUBJECTS_DIR,'fsaverage_sym','surf',[hemis{hh} '.' outmaps{m} '.avg.sym.' outName '.mgh']);
            else % prepare the data on the subject's native surface
                infile = fullfile(session_dir,[hemis{hh} '.cortex.avg.' inmaps{m} '.prfs.nii.gz']);
                outfile = fullfile(SUBJECTS_DIR,subject_name,'surf',[hemis{hh} '.' outmaps{m} '.avg.mgh']);
            end
            tmpfile = fullfile(session_dir,[hemis{hh} '.tmp.prfs.nii.gz']);
            % Convert radians to degrees (and Mathematica convention 0 - 180)
            if strcmp(inmaps{m},'copol')
                tmp = load_nifti(infile);
                % Flip like left hemisphere
                if strcmp(hemis{hh},'rh') && sy == 1 % sym space
                    upper = tmp.vol>0;
                    lower = tmp.vol<0;
                    tmp.vol(upper) = -(tmp.vol(upper) - pi);
                    tmp.vol(lower) = -(tmp.vol(lower) + pi);
                elseif strcmp(hemis{hh},'rh') && sy == 2 % native space
                    upper = tmp.vol>0;
                    lower = tmp.vol<0;
                    tmp.vol(upper) = tmp.vol(upper) - pi;
                    tmp.vol(lower) = tmp.vol(lower) + pi;
                end
                % Convert
                tmp.vol = rad2deg(tmp.vol) + 90;
                tmp.vol(tmp.vol<0) = nan;
                tmp.vol(tmp.vol>180) = nan;
                save_nifti(tmp,tmpfile);
                % Convert to mgh for Mathematica
                system(['mri_convert ' tmpfile ' ' outfile]);
                system(['rm ' tmpfile]);
            else
                % Convert to mgh for Mathematica
                system(['mri_convert ' infile ' ' outfile]);
                %             end
            end
        end
    end
end
disp('done.');