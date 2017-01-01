function copy_best_template(session_dir,template)

% copies the best fine schira model template ('ecc','pol','areas') to the
%   session_dir
%
%   usage:
%   copy_best_fine_template(session_dir)
%
%   Written by Andrew S Bock Oct 2015

%% set defaults
if ~exist('template','var')
    template = 'fine';
end
hemis = {'lh' 'rh'};
maps = {'ecc' 'pol' 'areas'};
%% Get the best template
for hh = 1:length(hemis)
    hemi = hemis{hh};
    [~,~,sorted_templates] = find_best_template(session_dir,template,hemi);
    % pull out the name
    varind = strfind(sorted_templates{1},'varexp');
    best_template = sorted_templates{1}(4:varind-2);
    for mm = 1:length(maps)
        switch template
            case 'coarse'
                invol = fullfile(session_dir,'pRFs','model_templates',...
                    [hemi '.' maps{mm} '.' best_template '.nii.gz']);
                outvol = fullfile(session_dir,...
                    [hemi '.' maps{mm} '.coarse.nii.gz']);
            case 'fine'
                invol = fullfile(session_dir,'pRFs','fine_model_templates',...
                    [hemi '.' maps{mm} '.' best_template '.nii.gz']);
                outvol = fullfile(session_dir,...
                    [hemi '.' maps{mm} '.fine.nii.gz']);
        end
        [~,~] = system(['cp ' invol ' ' outvol]);
    end
end