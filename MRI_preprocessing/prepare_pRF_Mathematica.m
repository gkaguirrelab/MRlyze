function prepare_pRF_Mathematica(session_dir,subject_name)

% Converts .nii.gz pRF cortical maps that result from 'average_pRF.m' to
% .mgh format, as well as changes the naming convention to that used by
% Mathematica pRF template fitting notebook
%
%   Usage:
%   prepare_pRF_Mathematica(session_dir,subject_name,sym)
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
                outfile = fullfile(SUBJECTS_DIR,subject_name,'surf',[hemis{hh} '.' outmaps{m} '.avg.sym.mgh']);
            else % prepare the data on the subject's native surface
                infile = fullfile(session_dir,[hemis{hh} '.cortex.avg.' inmaps{m} '.prfs.nii.gz']);
                outfile = fullfile(SUBJECTS_DIR,subject_name,'surf',[hemis{hh} '.' outmaps{m} '.avg.mgh']);
            end
            tmpfile = fullfile(session_dir,[hemis{hh} '.tmp.prfs.nii.gz']);
            % Convert radians to degrees (and Mathematica convention 0 - 180)
            if strcmp(inmaps{m},'copol')
                tmp = load_nifti(infile);
                % Flip like left hemisphere (but upside down)
                if strcmp(hemis{hh},'rh')
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