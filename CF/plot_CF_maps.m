function plot_CF_maps(session_dir,map_type,srcROI,trgROI,thresh,hemis)
%% Set defaults
if ~exist('srcROI','var')
    srcROI = 'volume';
end
if ~exist('trgROI','var')
    trgROI = 'prf_V1';
end
if ~exist('hemis','var')
    hemis = {'lh' 'rh'};
end
if ~exist('thresh','var')
    thresh = 0.125;
end
%% ROI mask
disp('Loading average maps...');
% Gaussian params
for s = 1:4 % assumes 4 Gaussian parameters
    if strcmp(srcROI,'volume');
        eval(['mh.sig' num2str(s) ' = load_nifti(fullfile(session_dir,''' ...
            'mh.' srcROI '.' trgROI '.' map_type '.avg.varsig' num2str(s) '.mni.nii.gz''' '));']);
    elseif strcmp(srcROI,'cortex')
        eval(['mh.sig' num2str(s) ' = load_nifti(fullfile(session_dir,''' ...
            'mh.' srcROI '.' trgROI '.' map_type '.avg.varsig' num2str(s) '.sym.nii.gz''' '));']);
    end
end
% variance explained
if strcmp(srcROI,'volume');
    eval(['mh.var = load_nifti(fullfile(session_dir,''' 'mh.' srcROI ...
        '.' trgROI '.' map_type '.avg.var.mni.nii.gz''' '));']);
elseif strcmp(srcROI,'cortex')
    eval(['mh.var = load_nifti(fullfile(session_dir,''' 'mh.' srcROI ...
        '.' trgROI '.' map_type '.avg.var.sym.nii.gz''' '));']);
end
for hh=1:length(hemis);
    hemi = hemis{hh};
    if strcmp(srcROI,'volume');
        eval([hemi '.areas = load_nifti(fullfile(session_dir,''' hemi '.LGN_SC.mni.' map_type '.nii.gz''' '));']);
    elseif strcmp(srcROI,'cortex')
        eval([hemi '.areas = load_nifti(fullfile(session_dir,''' 'lh.areas_pRF.sym.nii.gz''' '));']);
    end
    %% Find good voxels
    eval([hemi '.isGoodFit = mh.var.vol > ' num2str(thresh) ';']);
    if strcmp(srcROI,'volume')
        eval([hemi '.isLGN = ' hemi '.areas.vol == 1;']);
        eval([hemi '.isSC = ' hemi '.areas.vol == 2;']);
        eval([hemi '.goodLGN = ' hemi '.isLGN & ' hemi '.isGoodFit;']);
        eval([hemi '.goodSC = ' hemi '.isSC & ' hemi '.isGoodFit;']);
    elseif strcmp(srcROI,'cortex')
        eval([hemi '.isV2 = ' hemi '.areas.vol==2 | ' hemi '.areas.vol ==-2;']);
        eval([hemi '.isV3 = ' hemi '.areas.vol==3 | ' hemi '.areas.vol ==-3;']);
        eval([hemi '.goodV2 = ' hemi '.isV2 & ' hemi '.isGoodFit;']);
        eval([hemi '.goodV3 = ' hemi '.isV3 & ' hemi '.isGoodFit;']);
    end
end
disp('done.');
%% Get good values
for s = 1:4
    if strcmp(srcROI,'volume')
        % ipsi
        eval(['mh.sigvalsLGNI(:,' num2str(s) ') = mh.sig' num2str(s) ...
            '.vol(lh.goodLGN);']);
        eval(['mh.sigvalsSCI(:,' num2str(s) ') = mh.sig' num2str(s) ...
            '.vol(lh.goodSC);']);
        eval(['mh.sigvalsLGNC(:,' num2str(s) ') = mh.sig' num2str(s) ...
            '.vol(rh.goodLGN);']);
        eval(['mh.sigvalsSCC(:,' num2str(s) ') = mh.sig' num2str(s) ...
            '.vol(rh.goodSC);']);
        
    elseif strcmp(srcROI,'cortex')
        % V2
        eval(['mh.sigvalsV2(:,' num2str(s) ') = mh.sig' num2str(s)...
            '.vol(lh.goodV2);']);
        % V3
        eval(['mh.sigvalsV3(:,' num2str(s) ') = mh.sig' num2str(s)...
            '.vol(lh.goodV3);']);
    end
end
%% Save data file
save(fullfile(session_dir,['CF_data.' srcROI '.' trgROI '.' map_type '.mat']),'mh','lh','rh');
%% Plot Gaussians
if strcmp(srcROI,'volume')
    region = {'SCI' 'SCC' 'LGNI' 'LGNC'};
    regionTitle = {'Ipsi SC' 'Contra SC' 'Ipsi LGN' 'Contra LGN'};
elseif strcmp(srcROI,'cortex')
    region = {'V2' 'V3'};
    regionTitle = {'V2' 'V3'};
end
GaussianTitle = {'Center - sigma' 'Surround - sigma' 'Center - amplitude' 'Surround - amplitude'};
GaussianLims = {[0 10] [0 10] [-1 1] [-1 1]};
GaussianTics = {0:1:10 0:1:10 -1:.1:1 -1:.1:1};
for i = 1:length(region)
    if strcmp(srcROI,'volume') && strcmp(region{i},'LGNC')
        if ~exist(fullfile(session_dir,[hemi '.LGN.mni.' map_type '.nii.gz']),'file');
            disp('No rh.LGN.mni.nii.gz');
            return;
        end
    end
    ROI = [];
    eval(['ROI = [mh.sigvals' region{i} '];'])
    CovB=cov(ROI);
    beta = mean(ROI);
    R=ROI-(repmat(beta,size(ROI,1),1));
    MSE=mean(mean(R.^2)); % mse of R
    modelfun = @(b,x)(b(3)*(exp((-(x.^2))/(2*b(1).^2))) - b(4)*(exp((-(x.^2))/(2*b(2).^2))));
    xrange = (-25:.5:25)';
    [ypred,delta] = nlpredci(modelfun,xrange,beta,R,'Covar',CovB,...
        'MSE',MSE,'SimOpt','off');
    figure;h=plot(xrange,ypred,xrange,ypred-delta,xrange,ypred+delta);
    xlim([min(xrange) max(xrange)]);ylim([-1 1]);
    set(h(1),'linewidth',3,'Color','k','LineStyle','-');
    set(h(2),'linewidth',1,'Color','r','linestyle','--');
    set(h(3),'linewidth',1,'Color','r','linestyle','--');
    title([regionTitle{i} ' - V1'],'FontSize',30);
    axis square;
    set(gca,'TickDir','out')
    data = ROI;
    %mROI = mean(ROI);
    %         for d=1:4
    %             data(:,d) = data(:,d)/mROI(d);
    %         end
    %[coeff,score,latent,tsquared,explained,mu] = pca(data,'Centered',false);
    [coeff,score,latent,tsquared,explained,mu] = pca(data);
    %
    for pc=1%:2
        sig = (coeff(:,pc)+mu')';
        tmppc = make_DoG(xrange,sig);
        figure;h=plot(xrange,tmppc);
        xlim([min(xrange) max(xrange)]);ylim([-1 1]);
        title([regionTitle{i} ' - V1; PC' num2str(pc) ' - ' num2str(explained(pc)) ' variance explained'],'FontSize',10);
        xlabel('Distance (mm)','FontSize',20);
        ylabel('Amplitude','FontSize',20);
        axis square;
        set(gca,'TickDir','out')
    end
    %figure;[counts,centers] = hist(ROI,100);h=bar(centers,counts,'stacked','barwidth',1);axis square
    for s = 1:size(ROI,2);
        figure;hist(ROI(:,s),50);
        xlim(GaussianLims{s});set(gca,'XTick',GaussianTics{s});
        title([regionTitle{i} ' - V1; ' GaussianTitle{s}],'FontSize',20);
        set(gca,'TickDir','out')
        xlabel('Gaussian parameter value','FontSize',20);
        ylabel('Number of Voxels/Vertices','FontSize',20);
    end
end