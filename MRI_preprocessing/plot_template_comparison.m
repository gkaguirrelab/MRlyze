function mat = plot_template_comparison(session_dirs,templates,manual_templates,hemis)

% Makes a bar plot, comparing variance explained by different template
% types.
%
%   Written by Andrew S Bock Oct 2015

%% set defaults
if ~exist('session_dirs','var')
    session_dirs = {...
        '/data/jet/abock/data/Retinotopy/AEK/10012014' ...
        '/data/jet/abock/data/Retinotopy/ASB/10272014' ...
        '/data/jet/abock/data/Retinotopy/GKA/10152014' ...
        };
end
if ~exist('templates','var')
    templates = {...
        'pRF' ...
        'fine' ...
        'coarse' ...
        'manual' ...
        'anat' ...
        };
end
if ~exist('hemis','var')
    hemis = {'lh' 'rh'};
end
% manual_templates = {...
%     '0.6,-0.2,-0.5' ...
%     '0.0,-0.3,0.2' ...
%     '0.4,-0.3,-0.3' ...
%     '0.0,-0.2,0.2' ...
%     '0.4,0.1,-0.4' ...
%     '0.3,-0.1,-0.1' ...
%     };
if ~exist('manual_templates','var')
    % AEK
    manual_templates{1,1} = '10.6.4';
    manual_templates{1,2} = '4.5.11';
    % ASB
    manual_templates{2,1} = '8.5.6';
    manual_templates{2,2} = '4.6.11';
    % GKA
    manual_templates{3,1} = '8.9.5';
    manual_templates{3,2} = '7.7.8';
end
%% Pull out the template values
mat = nan(length(templates),length(session_dirs),length(hemis)); % matrix
for tt = 1:length(templates)
    template = templates{tt};
    for ss = 1:length(session_dirs)
        session_dir = session_dirs{ss};
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            switch template
                case 'pRF'
                    if strcmp(session_dir,'/data/jet/abock/data/Retinotopy/AEK/10012014')
                        runs = [3,4,6];
                    else
                        runs = [2,4,6];
                    end
                    tmp = regress_template(session_dir,runs,hemi,template,0);
                    mat(tt,ss,hh) = nansum(tmp(2:3));
                case 'fine'
                    tmp = find_best_template(session_dir,template,hemi);
                    mat(tt,ss,hh) = tmp(1);
                case 'coarse'
                    tmp = find_best_template(session_dir,template,hemi);
                    mat(tt,ss,hh) = tmp(1);
                case 'manual'
                    if strcmp(session_dir,'/data/jet/abock/data/Retinotopy/AEK/10012014')
                        runs = [3,4,6];
                    else
                        runs = [2,4,6];
                    end
                    man_temp = manual_templates{ss,hh};
                    tmp = regress_template(session_dir,runs,hemi,man_temp,0);
                    mat(tt,ss,hh) = nansum(tmp(2:3));
                case 'anat'
                    if strcmp(session_dir,'/data/jet/abock/data/Retinotopy/AEK/10012014')
                        runs = [3,4,6];
                    else
                        runs = [2,4,5];
                    end
                    tmp = regress_template(session_dir,runs,hemi,template,0);
                    mat(tt,ss,hh) = nansum(tmp(2:3));
            end
        end
    end
end
%% Make a bar plot
newmat = (reshape(mat,size(mat,1),size(mat,2)*size(mat,3)))';
sem = std(newmat) / sqrt(size(newmat,1));
figure;bar(mean(newmat));hold on;
errorbar(mean(newmat),sem,'.','MarkerSize',0.0001);
ylim([0 0.05]);
xlabel('Template Type','FontSize',20);
ylabel('Variance Explained','FontSize',20);
ax = gca;
set(ax,'XTickLabel',templates,'FontSize',15);
axis square;