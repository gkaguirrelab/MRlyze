function save_template_correlation_matrices(session_dir,hemi,func,cond,runs,templateSize,V2V3)

% Saves both the observed and predicted the cross-correlation matrices
% which result from the retinotopic template fitting pipeline.
%
%   Written by Andrew S Bock Jan 2016

%% Set defaults
if ~exist('hemi','var')
    hemi = 'lh';
end
if ~exist('func','var')
    func = 'wdrf.tf';
end
if ~exist('cond','var')
    cond = 'Movie';
end
if ~exist('runs','var')
    runs = [2 4 6];
end
if ~exist('templateSize','var')
    templateSize = 'full';
end
if V2V3
    saveName = 'V2V3';
else
    saveName = 'V1';
end
%% Define variables
% templates{1} = {...
%     'pRF'};
% templates{1} = {'4.6.4' '5.2.7'}; % last two are 'best' and 'worst', repsectively for GKA lh
% templateTypes = {'coarse'};
templateTypes = {'pRF' 'coarse'};
templates{1} = {'pRF'}; 
templates{2} = {'4.6.4' '5.2.7'}; % 'best' and 'worst', repsectively for GKA lh
%templateTypes = {'pRF' 'coarse' 'fine'};
% Create custome white -> blue color scale
blueScale = zeros(100,3);
blueScale(:,1) = linspace(1,0.4,100); % transition to white -> blue
blueScale(:,2) = linspace(1,0.4,100); % transition to white -> blue
blueScale(:,3) = ones(100,1); % transition to white -> blue
% Create custome white -> red color scale
redScale = zeros(100,3);
redScale(:,1) = ones(100,1); % set red
redScale(:,2) = linspace(1,0.4,100); % transition to white -> red
redScale(:,3) = linspace(1,0.4,100); % transition to white -> red
% Create custome white -> gray color scale
grayScale = zeros(100,3);
grayScale(:,1) = linspace(1,0.4,100); % transition to white -> gray
grayScale(:,2) = linspace(1,0.4,100); % transition to white -> gray
grayScale(:,3) = linspace(1,0.4,100); % transition to white -> gray
%% Loop through templates, save images
for tt = 1:length(templateTypes)
    templateType = templateTypes{tt};
    for t = 1:length(templates{tt})
        template = templates{tt}{t};
        switch templateType
            case 'pRF'
                pRF_dir = fullfile(session_dir,'pRFs','pRF_templates','decimated_templates');
            case 'coarse'
                pRF_dir = fullfile(session_dir,'pRFs','coarse_model_templates','decimated_templates');
            case 'fine'
                if V2V3
                    pRF_dir = fullfile(session_dir,'pRFs','fine_model_templates','V2V3','decimated_templates');
                else
                    pRF_dir = fullfile(session_dir,'pRFs','fine_model_templates','V1','decimated_templates');
                end
        end
        % Get variance explained
        if V2V3
            tdir = fullfile(session_dir,'pRFs',templateType,func,cond,'V2V3');
            varexp = load(fullfile(tdir,[hemi '.' template '.varexp.txt']));
            varexp = nansum(varexp(1:3)); %%% V1-V2, V1-V3, AND V2-V3 %%%%
        else
            tdir = fullfile(session_dir,'pRFs',templateType,func,cond,'V1');
            varexp = load(fullfile(tdir,[hemi '.' template '.varexp.txt']));
            varexp = nansum(varexp(1:2)); %%% V1-V2, V1-V3 %%%%
        end
        % Create the observed and predicted matrices
        [obs_all,pred_all,~,~,rowV1inds,rowV2inds,rowV3inds,colV1inds,colV2inds,colV3inds] = create_template_matrices(...
            session_dir,pRF_dir,template,hemi,func,runs,templateSize);
        %% Show Movie
        figure('units','normalized','position',[0 0 1 1]);
        imagesc(fliplr(obs_all'),[0.25 1]);
        colormap(blueScale);
        title(['Observed ''' template ''' Correlation Matrix'],'FontSize',30);
        xlabel('Visual Areas','FontSize',20);
        switch templateSize
            case 'full'
                set(gca,'XTick',[...
                    size(obs_all,1) - median(rowV3inds),...
                    size(obs_all,1) - median(rowV2inds),...
                    size(obs_all,1) - median(rowV1inds)]);
                set(gca,'XTickLabel',{'V3','V2','V1'},'FontSize',15);
                axis square
            case 'V2_V3'
                if V2V2
                    set(gca,'XTick',[size(obs_all,1) - median(rowV2inds),...
                        size(obs_all,1) - median(rowV1inds)]);
                    set(gca,'XTickLabel',{'V2','V1'},'FontSize',15);
                else
                    set(gca,'XTick',median(rowV1inds));
                    set(gca,'XTickLabel',{'V1'},'FontSize',15);
                end
        end
        ylabel('Visual Areas','FontSize',20);
        switch templateSize
            case 'full'
                set(gca,'YTick',[median(colV1inds),median(colV2inds),median(colV3inds)]);
                set(gca,'YTickLabel',{'V1','V2','V3'},'FontSize',15);
            case 'V2_V3'
                set(gca,'YTick',[median(colV2inds),median(colV3inds)]);
                set(gca,'YTickLabel',{'V2','V3'},'FontSize',15);
        end
        hcb=colorbar;
        hL = ylabel(hcb,'Correlation (fisher-z)','FontSize',20);
        set(hL,'Rotation',90);
        savefigs('pdf',[templateType '_' template '_observed_corr_mat_' templateSize '_' saveName]);
        close all;
        %% Show Template Prediction
        figure('units','normalized','position',[0 0 1 1]);
        imagesc(fliplr(pred_all'),[0.25 1]);
        colormap(grayScale);
        title(['Predicted ''' template ''' Correlation Matrix - Variance explained = ' ...
            num2str(varexp)],'FontSize',30);
        xlabel('Visual Areas','FontSize',20);
        switch templateSize
            case 'full'
                set(gca,'XTick',[...
                    size(obs_all,1) - median(rowV3inds),...
                    size(obs_all,1) - median(rowV2inds),...
                    size(obs_all,1) - median(rowV1inds)]);
                set(gca,'XTickLabel',{'V3','V2','V1'},'FontSize',15);
                axis square
            case 'V2_V3'
                if V2V3
                    set(gca,'XTick',[size(obs_all,1) - median(rowV2inds),...
                        size(obs_all,1) - median(rowV1inds)]);
                    set(gca,'XTickLabel',{'V2','V1'},'FontSize',15);
                else
                    set(gca,'XTick',median(rowV1inds));
                    set(gca,'XTickLabel',{'V1'},'FontSize',15);
                end
        end
        ylabel('Visual Areas','FontSize',20);
        switch templateSize
            case 'full'
                set(gca,'YTick',[median(colV1inds),median(colV2inds),median(colV3inds)]);
                set(gca,'YTickLabel',{'V1','V2','V3'},'FontSize',15);
            case 'V2_V3'
                set(gca,'YTick',[median(colV2inds),median(colV3inds)]);
                set(gca,'YTickLabel',{'V2','V3'},'FontSize',15);
        end
        hcb=colorbar;
        hL = ylabel(hcb,'Correlation (fisher-z)','FontSize',20);
        set(hL,'Rotation',90);
        savefigs('pdf',[templateType '_' template '_predicted_corr_mat_' templateSize '_' saveName]);
        close all;
    end
end