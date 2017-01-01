function plot_pRF_surf(session_dir,subject,func,roi)

% show pRF and ccpRF estimates on surface (patch); if simulated, also show estimated vs. true
% PB 03/2013
% modified AB 03/2013

%% Set defaults
if ~exist('session_dir','var')
    error('No ''session_dir'' defined')
end
if ~exist('subject','var')
    error('No ''subject'' defined')
end
if ~exist('func','var')
    func = {'sdbrf'}; % prefix of functional files % note the .tf was left off for Ari's pRF analysis files
end

if ~exist('roi','var')
    roi = 3;% ROI 1=V1;2=V2;3=occipital;
end
analysis_type = 'pRF';
nruns = 3;
%% Set up directories and variables
prf_dir = fullfile(session_dir,'prfs');
d = listdir(fullfile(prf_dir,'pRF*'),'dirs');
prf_dir = fullfile(prf_dir,d{1});
save_dir = session_dir;
hemi = {'lh' 'rh'};
rois = {'V1' 'V1_V3' 'occipital'};%{'V1' 'V2' 'V3' 'brain'}; % for ccpRF analysis, rois{1} must be V1
% inclusion crtieria
maps = {'ecc'; 'polang'; 'sig'; 'co'};%; 'amp' ;'mse'
thresh{1} = [-inf inf];
thresh{2} = [-inf inf];
thresh{3} = [-inf inf];
thresh{4} = [-inf inf];
%% Convert values to myfields
for dd = 1:length(func);
    for h = 1:length(hemi)
        for r = 1:nruns
            % find the appropriate fit mat file
            cd(fullfile(prf_dir,[func{dd} '.' hemi{h} '.' num2str(r)]));
            matfile = listdir('./*fFit*','files');
            % Load in mat file, extract parameters
            F = load(matfile{1});
            model = F.model;
            fit = 1 - (model{1}.rss ./ model{1}.rawrss);
            fit(~isfinite(fit)) = 0;
            fit = max(fit, 0);fit = min(fit, 1);
            fit(fit == 0) = NaN;
            [theta, rho] = cart2pol(model{1}.x0, model{1}.y0);
            % Load up a file in the same subject surface space
            nii = load_nifti(fullfile(session_dir,[hemi{h} '.ecc.nii.gz']));
            % Save nifti
            % note the .tf was left off for Ari's pRF analysis files
            for m = 1:length(maps)
                savename=fullfile(save_dir,[hemi{h} '_' func{dd} '.tf_' ...
                    analysis_type '_' rois{roi} '_run' num2str(r) '_' maps{m}]);
                switch maps{m}
                    case 'ecc'
                        nii.vol = rho';
                        save_nifti(nii,[savename '_surf.nii.gz']);
                    case 'polang'
                        nii.vol = -theta';
                        if h == 2
                            pos=find(nii.vol>0);
                            neg=find(nii.vol<0);
                            nii.vol(pos)=-nii.vol(pos) + pi;
                            nii.vol(neg)=-nii.vol(neg) - pi;
                        end
                        save_nifti(nii,[savename '_surf.nii.gz']);
                    case 'sig'
                        nii.vol = model{1}.sigma.major';
                        save_nifti(nii,[savename '_surf.nii.gz']);
                    case 'co'
                        nii.vol = fit';
                        save_nifti(nii,[savename '_surf.nii.gz']);
                end
                if h ==1
                    [~,~] = system(['mri_surf2surf --hemi ' hemi{h} ' --srcsubject ' subject ' --srcsurfval ' savename '_surf.nii.gz --trgsubject fsaverage_sym --trgsurfval ' savename '_fsavgsurf.nii.gz']);
                else
                    [~,~] = system(['mri_surf2surf --hemi lh --srcsubject ' subject '/xhemi --srcsurfval ' savename '_surf.nii.gz --trgsubject fsaverage_sym --trgsurfval ' savename '_fsavgsurf.nii.gz']);
                end
            end
        end
    end
end