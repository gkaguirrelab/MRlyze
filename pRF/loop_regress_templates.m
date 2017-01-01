function loop_regress_templates(session_dir,runs,hemi,templateType,cluster)

% runs a for loop through the templates in the specified session_dir and
% template type ('coarse' or 'fine'), running either regress_template or
% regress_fine_template, respectively.
%
%   Usage:
%   loop_regress_templates(session_dir,runs,hemi,templateType)
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
if ~exist('cluster','var')
    cluster = 1; % assume running on the cluster
end
%%
switch templateType
    case 'coarse'
        templateDir = fullfile(session_dir,'pRFs','model_templates','decimated_templates');
        templates = listdir(fullfile(templateDir,[hemi '.areas.*.nii.gz']),'files');
        progBar = ProgressBar(length(templates),'Looping through templates...');
        for i = 1:length(templates)
            % pull out the template name
            niidot = strfind(templates{i},'.nii.gz');
            template = templates{i}(10:niidot-1);
            % run regress_template
            regress_template(session_dir,runs,hemi,template,cluster);
            progBar(i);
        end
    case 'fine'
        templateDir = fullfile(session_dir,'pRFs','fine_model_templates','decimated_templates');
        templates = listdir(fullfile(templateDir,[hemi '.areas.*.nii.gz']),'files');
        progBar = ProgressBar(length(templates),'Looping through templates...');
        for i = 1:length(templates)
            % pull out the template name
            niidot = strfind(templates{i},'.nii.gz');
            template = templates{i}(10:niidot-1);
            % run regress_fine_template
            regress_fine_template(session_dir,runs,hemi,template,cluster);
            progBar(i);
        end
end