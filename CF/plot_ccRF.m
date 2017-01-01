function plot_ccRF(session_dir,subject_name,rr,func,template,roi,hemi)

% plot pRF and ccpRF estimates on surface (patch); if simulated, also show estimated vs. true
% PB 03/2013
% modified AB 03/2013

if ~exist('session_dir','var')
    error('No ''session_dir'' defined')
end
if ~exist('subject_name','var')
    error('No ''subject_name'' defined')
end
if ~exist('rr','var')
    error('Run ''rr'' variable not defined')
end
if ~exist('func','var')
    func = 'sdbrf.tf'; % prefix of functional files
end
if ~exist('template','var')
    template = 'anatomical';
end
if ~exist('roi','var')
    roi = 3;% occipital
end
if ~exist('hemi','var')
    hemi = {'lh' 'rh'};
end
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
d = listdir(fullfile(session_dir,'*bold_*'),'dirs');
nruns = length(d);
save_dir = session_dir;
rois = {'V1' 'V1_V3' 'occipital' 'cortex' 'subcortical'};% for ccpRF analysis, rois{1} must be V1
maps = {'x' 'y' 'ecc' 'pol' 'sig' 'co'};
templates = {'ecc' 'pol' 'areas'};
% inclusion crtieria
myfields = {'x'; 'y'; 'ecc'; 'polang'; 'sig'; 'co'};%; 'amp' ;'mse'
for iv = 1:length(myfields)
    mythresholds{iv,1} = [-Inf +Inf];
    plotlims{iv,1} = [-10 10];
end
% mythresholds{3} = [0 90];
% mythresholds{4} = [-pi pi];
% mythresholds{5} = [0 inf];
% mythresholds{6} = [0 1];

mythresholds{3} = [-inf inf];
mythresholds{4} = [-inf inf];
mythresholds{5} = [-inf inf];
mythresholds{6} = [-inf inf];

plotlims{3} = [0 16];
plotlims{4} = [-pi pi];
plotlims{5} = [0 15];
plotlims{6} = [mythresholds{6}(1)*0 1];
params = [myfields,mythresholds,plotlims];
%% Plot in the volume
if roi == 5
    goodprfs = cell(2,1);
    disp('Plotting pRFs in volume...');
    system(['fslroi ' fullfile(session_dir,d{rr},[func '.nii.gz']) ' ' ...
        fullfile(session_dir,d{rr},'tmp.nii.gz') ' 0 1']);
    for hh=1:length(hemi) %hemispheres
        ROI_file = fullfile(session_dir,d{rr},'tmp.nii.gz'); % load up a volume nii to be used as ROI nii
        % load prfs & brain surface
        input_filename = fullfile(session_dir,'prfs',[func '_' template],[hemi{hh} '_' func ...
            '_' template '_' rois{roi} '_run' num2str(rr)  '.mat']);
        load(input_filename)
        % store anatomical vars
        if strcmp(template,'Stim') % if stimulus drive and do_pRF
            ver4plot{hh} = vc2vx.voxel;
            N{hh} = length(vc2vx.voxel);
        else
            ver4plot{hh} = Vertices4plot{roi,hh}.voxel;
            N{hh} = length(Vertices4plot{roi,hh}.voxel);
        end
        % rearrange prfs struct
        tmp = cat(1,pRFs.center);
        p.x = tmp(:,1);
        p.y = tmp(:,2);
        switch template
            case 'prf'
                p.x = s.x(tmp(:,1));
                p.y = s.y(tmp(:,1));
                p.ecc = s.ecc(tmp(:,1));
                p.polang = s.pol(tmp(:,1));
            case 'anatomical'
                p.x = s.x(tmp(:,1));
                p.y = s.y(tmp(:,1));
                p.ecc = s.ecc(tmp(:,1));
                p.polang = s.pol(tmp(:,1));
            case 'Templ'
                [p.polang, p.ecc] = cart2pol(p.x,p.y);
            case 'Stim'
                [p.polang, p.ecc] = cart2pol(p.x,p.y);
        end
        for iv = 1:length(myfields)
            if isfield(pRFs,params{iv,1})
                eval(['p.' params{iv,1} ' = cat(1,pRFs.' params{iv,1} ');'])
                eval(['p.' params{iv,1} ' = p.' params{iv,1} '(:,1);'])
            end
        end
        % thresholding
        failed = zeros(N{hh},length(myfields));
        for iv = 1:size(params,1)
            eval(['thisvar = p.' params{iv,1} ';']);
            failed(:,iv) = thisvar < params{iv,2}(1) | thisvar > params{iv,2}(2) | isnan(thisvar);
        end
        % index to prfs that pass all the thresholds
        goodprfs{hh} = sum(failed,2) == 0;
        
        for iv = 1:size(params,1)
            % set to nan all non selected
            eval(['p.' params{iv,1} '(~goodprfs{hh}) = NaN;']);
            % also set to nan all non selected in the first condition (fullfieROI_surf_fileld)
            eval(['p.' params{iv,1} '(~goodprfs{hh}) = NaN;']);
        end
        % store prfs
        prfs{hh} = p;
        % Plot maps on native and fsaverage_sym surface
        for iv = 3:length(myfields) %ecc, polang, sig, co
            thisvar = [];
            ind = thisvar;
            vals = thisvar;
            eval(['thisvar = prfs{hh}.' params{iv,1} ';'])
            
            % make surface maps
            for i = 1:N{hh}
                % select all ver4plot that project to this voxel
                vc4thisvox = ver4plot{hh}(i).vertices;
                ind = [ind; vc4thisvox];
                vals = [vals; ones(length(vc4thisvox),1)*thisvar(i)];
            end
            %create nifti file for each map (ecc, pol, sig, co)
            newROInii=load_nifti(ROI_file);
            newimg=nan(size(newROInii.vol));
            newimg(ind)= vals;
            newROInii.vol=newimg; % replace the image
            savename=fullfile(save_dir,[hemi{hh} '_' func '_' ...
                template '_' rois{roi} '_run' num2str(rr) ...
                '_' maps{iv}]);
            save_nifti(newROInii,[savename '_vol.nii.gz']);
            [~,~] = system(['mri_vol2vol --mov ' savename '_vol.nii.gz ' ...
                '--targ ' fullfile(SUBJECTS_DIR,subject_name,'mri','orig.mgz') ...
                ' --o ' savename '_anatvol.nii.gz --reg ' ...
                fullfile(session_dir,d{rr},'brf_bbreg.dat')]);
        end
    end
end
%% Threshold and transform to fsaverage surface
if roi ~=5
    goodprfs = cell(2,1);
    disp('Plotting pRFs on surface...');
    for hh=1:length(hemi) %hemispheres
        ROI_surf_file = fullfile(session_dir,[hemi{hh} '.ecc.nii.gz']); % load up a volume nii to be used as ROI nii
        % load prfs & brain surface
        
        input_filename = fullfile(session_dir,'prfs',[func '_' template],[hemi{hh} '_' func ...
            '_' template '_' rois{roi} '_run' num2str(rr)  '.mat']);
        
        load(input_filename)
        
        % store anatomical vars
        if strcmp(template,'Stim') % if stimulus drive and do_pRF
            ver4plot{hh} = vc2vx.voxel;
            N{hh} = length(vc2vx.voxel);
        else
            ver4plot{hh} = Vertices4plot{roi,hh}.voxel;
            N{hh} = length(Vertices4plot{roi,hh}.voxel);
        end
        % rearrange prfs struct
        tmp = cat(1,pRFs.center);
        p.x = tmp(:,1);
        p.y = tmp(:,2);
        
        switch template
            case 'prf'
                p.x = s.x(tmp(:,1));
                p.y = s.y(tmp(:,1));
                p.ecc = s.ecc(tmp(:,1));
                p.polang = s.pol(tmp(:,1));
            case 'anatomical'
                p.x = s.x(tmp(:,1));
                p.y = s.y(tmp(:,1));
                p.ecc = s.ecc(tmp(:,1));
                p.polang = s.pol(tmp(:,1));
            case 'Templ'
                [p.polang, p.ecc] = cart2pol(p.x,p.y);
            case 'Stim'
                [p.polang, p.ecc] = cart2pol(p.x,p.y);
        end
        
        for iv = 1:length(myfields)
            if isfield(pRFs,params{iv,1})
                eval(['p.' params{iv,1} ' = cat(1,pRFs.' params{iv,1} ');'])
                eval(['p.' params{iv,1} ' = p.' params{iv,1} '(:,1);'])
            end
        end
        
        % thresholding
        failed = zeros(N{hh},length(myfields));
        for iv = 1:size(params,1)
            eval(['thisvar = p.' params{iv,1} ';']);
            failed(:,iv) = thisvar < params{iv,2}(1) | thisvar > params{iv,2}(2) | isnan(thisvar);
        end
        % index to prfs that pass all the thresholds
        goodprfs{hh} = sum(failed,2) == 0;
        
        for iv = 1:size(params,1)
            % set to nan all non selected
            eval(['p.' params{iv,1} '(~goodprfs{hh}) = NaN;']);
            % also set to nan all non selected in the first condition (fullfieROI_surf_fileld)
            eval(['p.' params{iv,1} '(~goodprfs{hh}) = NaN;']);
        end
        % store prfs
        prfs{hh} = p;
        
        % Plot subject pRF templates on fsaverage_sym surface
        for t = 1:length(templates) % ecc, pol, sig
            if hh ==1
                [~,~] = system(['mri_surf2surf --hemi ' hemi{hh} ' --srcsubject ' ...
                    subject_name ' --srcsurfval ' fullfile(session_dir,[hemi{hh} '.' templates{t} '_pRF.nii.gz']) ...
                    ' --trgsubject fsaverage_sym --trgsurfval ' ...
                    fullfile(session_dir,[hemi{hh} '.' templates{t} '_pRF_fsavgsurf.nii.gz'])]);
            else
                [~,~] = system(['mri_surf2surf --hemi lh --srcsubject ' ...
                    subject_name '/xhemi --srcsurfval ' fullfile(session_dir,[hemi{hh} '.' templates{t} '_pRF.nii.gz']) ...
                    ' --trgsubject fsaverage_sym --trgsurfval ' ...
                    fullfile(session_dir,[hemi{hh} '.' templates{t} '_pRF_fsavgsurf.nii.gz'])]);
            end
        end
        % Plot maps on native and fsaverage_sym surface
        for iv = 3:length(myfields) %ecc, polang, sig, co
            thisvar = [];
            ind = thisvar;
            vals = thisvar;
            eval(['thisvar = prfs{hh}.' params{iv,1} ';'])
            
            % make surface maps
            for i = 1:N{hh}
                % select all ver4plot that project to this voxel
                vc4thisvox = ver4plot{hh}(i).vertices;
                ind = [ind; vc4thisvox];
                vals = [vals; ones(length(vc4thisvox),1)*thisvar(i)];
            end
            %create nifti file for each map (ecc, pol, sig, co)
            newROInii=load_nifti(ROI_surf_file);
            newimg=nan(size(newROInii.vol));
            newimg(ind)= vals;
            newROInii.vol=newimg; % replace the image
            if strcmp(maps(iv),'pol')
                newROInii.vol = newROInii.vol - pi/2;
            end
            savename=fullfile(save_dir,[hemi{hh} '_' func '_' ...
                template '_' rois{roi} '_run' num2str(rr) ...
                '_' maps{iv}]);            
            save_nifti(newROInii,[savename '_surf.nii.gz']);
            if hh ==1
                [~,~] = system(['mri_surf2surf --hemi ' hemi{hh} ' --srcsubject ' ...
                    subject_name ' --srcsurfval ' savename '_surf.nii.gz --trgsubject fsaverage_sym --trgsurfval ' ...
                    savename '_fsavgsurf.nii.gz']);
            else
                [~,~] = system(['mri_surf2surf --hemi lh --srcsubject ' ...
                    subject_name '/xhemi --srcsurfval ' savename '_surf.nii.gz --trgsubject fsaverage_sym --trgsurfval ' ...
                    savename '_fsavgsurf.nii.gz']);
            end
        end
    end
end