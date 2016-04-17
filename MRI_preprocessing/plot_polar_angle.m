function plot_polar_angle(session_dir)

% Plot polar angle by area in the Benson Template
%
%   usage:
%   plot_polar_angle(session_dir)

%% Set up defaults
if ~exist('session_dir','var')
    error('No ''session_dir'' defined')
end
tmp.areas = '~/data/2014-10-29.areas-template.nii.gz';
tmp.ecc = '~/data/2014-10-29.eccen-template.nii.gz';
tmp.polang = '~/data/2014-10-29.angle-template-RADS.nii.gz';
%% Load up template and pRF/ccRF
areas = load_nifti(tmp.areas);
ecc = load_nifti(tmp.ecc);
polang = load_nifti(tmp.polang);
polang.vol = rad2deg(polang.vol) + 90;
ccRFmovie = load_nifti(fullfile(session_dir,'mh_sdbrf.tf_ccRF_occipital_avg_movie_polang_fsavgsurf.nii.gz'));
ccRFmovie.vol = rad2deg(ccRFmovie.vol) + 90;
ccRFmovieco = load_nifti(fullfile(session_dir,'mh_sdbrf.tf_ccRF_occipital_avg_movie_co_fsavgsurf.nii.gz'));
ccRFbars = load_nifti(fullfile(session_dir,'mh_sdbrf.tf_ccRF_occipital_avg_bars_polang_fsavgsurf.nii.gz'));
pRFbars = load_nifti(fullfile(session_dir,'mh_sdbrf.tf_pRF_occipital_avg_polang_fsavgsurf.nii.gz'));
pRFbars.vol = rad2deg(pRFbars.vol) + 90;
pRFbarsco = load_nifti(fullfile(session_dir,'mh_sdbrf.tf_pRF_occipital_avg_co_fsavgsurf.nii.gz'));
%% Threshold by variance
var.ccRFmovie = (ccRFmovieco.vol).^2;
var.ccRFind = var.ccRFmovie<0.15;
var.pRFbars = pRFbarsco.vol;
var.pRFind = var.pRFbars<0.15;
ccRFmovie.vol(var.ccRFind) = nan;
pRFbars.vol(var.pRFind) = nan;
%% Find various indicies
ROI.V3v = areas.vol == 3;
ROI.V2v = areas.vol == 2;
ROI.V1v = areas.vol == 1;
ROI.V1d = areas.vol == -1;
ROI.V2d = areas.vol == -2;
ROI.V3d = areas.vol == -3;
% eccen.one = ecc.vol < 1;
% eccen.one_two = ecc.vol > 1 & ecc.vol < 2;
% eccen.two_four = ecc.vol > 2 & ecc.vol < 4;
% eccen.four_eight = ecc.vol > 4 & ecc.vol < 8;
% eccen.eight_sixteen = ecc.vol > 8 & ecc.vol < 16;
% eccen.sixteen_thirtytwo = ecc.vol > 16 & ecc.vol < 32;
%eccen.one_eight = ecc.vol > 1 & ecc.vol < 8;
eccen.two_eight = ecc.vol > 2 & ecc.vol < 8;

%% Plot polar angle by areas and eccentricty
Enames = fieldnames(eccen);
Rnames = fieldnames(ROI);
for E = 1:length(Enames);
    tvals = [];
    cvals = [];
    pvals = [];
    Xvals = 0;
    for R = 1:length(Rnames);
        ind = ROI.(Rnames{R}) & eccen.(Enames{E});
        tmpl = polang.vol(ind);
        ccRF = ccRFmovie.vol(ind);
        pRF = pRFbars.vol(ind);
        [B,idx] = sort(tmpl);
        tval = tmpl(idx);
        cval = ccRF(idx);
        pval = pRF(idx);
        if R == 1 || R == 3 || R == 4 || R == 6
            %    if R == 2 || R == 5
            tval = flipud(tval);
            cval = flipud(cval);
            pval = flipud(pval);
        end
        % Remove V1 from ccRF
        if R == 3 || R == 4
            cval = nan(size(cval));
        end
        Xvals = [Xvals;(length(idx)+Xvals(end))];
        tvals = [tvals;tval];
        cvals = [cvals;cval];
        pvals = [pvals;pval];
    end
    figure;
    %plot(tvals,'ko');hold on % plot the template curve
    % Add grey fill for every other visual region
%     for R = 1:length(Rnames)
%         if ~mod(R,2)
%             x1 = Xvals(R):Xvals(R+1);
%             y1 = zeros(size(x1));
%             y2 = 180*ones(size(x1));
%             X = [x1,fliplr(x1)];
%             Y = [y1,fliplr(y2)];
%             fill(X,Y,[.75 .75 .75]);hold on
%         end
%     end
    c = plot(cvals,'o'); hold on
    p = plot(pvals,'o'); hold on
    set(c,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5)
    set(p,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5)
    ylim([0 180]);
    xlim([0 Xvals(end)]);
end
%% Set Title, etc
T = title('Polar Angle: 2 - 8 degrees Ecc');
xl = xlabel('Visual Areas - V3v to V3d');
yl = ylabel('Polar Angle (radians)');
set(T,'FontSize',25)
set(xl,'Interpreter','none','FontSize',20);
set(yl,'Interpreter','none','FontSize',20);
% Remove tick marks
set(gca,'xtick',[])
set(gca,'ytick',[])
%% Plot eccentricty band
% ecc.vol(ecc.vol<2 | ecc.vol>8) = nan;
% save_nifti(ecc,'~/data/ecc_2_8.nii.gz')
% surface_plot('ecc','~/data/ecc_2_8.nii.gz','fsaverage_sym','lh')