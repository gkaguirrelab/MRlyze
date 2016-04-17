subjects = {'AEK' 'ASB' 'GKA'};
session_dirs = {...
    '/jet/abock/data/Retinotopy/AEK/10012014' ...
    '/jet/abock/data/Retinotopy/ASB/10272014' ...
    '/jet/abock/data/Retinotopy/GKA/06052015'};
ROIs = {'V2' 'V3' 'LGN' 'SC'};
hemis = {'lh' 'rh' ,'mh'};

for ss=1:length(subjects)
    for rr = 1:length(ROIs)
        for hh = 1:length(hemis)
            if strcmp(ROIs{rr},'LGN') || strcmp(ROIs{rr},'SC')
                tmp = get_HRF_peak(...
                    session_dirs{ss},hemis{hh},'volume','prf_V1','movie',ROIs{rr});
                datamat(ss,rr,hh,1) = tmp.cfmask;
                datamat(ss,rr,hh,2) = tmp.prfmask;
                datamat(ss,rr,hh,3) = tmp.cfmask_shift;
            else
                tmp = get_HRF_peak(...
                    session_dirs{ss},hemis{hh},'cortex','prf_V1','movie',ROIs{rr});
                datamat(ss,rr,hh,1) = tmp.cfmask;
                datamat(ss,rr,hh,2) = tmp.prfmask;
                datamat(ss,rr,hh,3) = tmp.cfmask_shift;
            end
        end
    end
end
%%
xlabels = {
    'V2 - HRF difference'; 'V2 - CF temporal offset'; ...
    'V3 - HRF difference'; 'V3 - CF temporal offset';...
    'LGN - HRF difference'; 'LGN - CF temporal offset';...
    'SC - HRF difference'; 'SC - CF temporal offset'};
for hh = 1:length(hemis)
    figure;
    ct = 0;
    for i = 1:4
        ct = ct + 1;
        data(:,ct) = datamat(:,i,hh,1) - datamat(:,i,hh,2);
        ct = ct + 1;
        data(:,ct) = datamat(:,i,hh,3);
    end
    boxplot(data,'labels',xlabels)
    xlabel('ROI','FontSize',30);
    ylabel('Temporal difference (seconds)','FontSize',30);
end
%%
rownames = {'AEK';'ASB';'GKA'};
xlabels = {
    'V2 - HRF in V1 CF'; 'V2 - HRF'; ...
    'V2 - HRF difference'; 'V2 - CF temporal offset'; ...
    'V3 - HRF in V1 CF'; 'V3 - HRF'; ...
    'V3 - HRF difference'; 'V3 - CF temporal offset';...
    'LGN - HRF in V1 CF'; 'LGN - HRF'; ...
    'LGN - HRF difference'; 'LGN - CF temporal offset';...
    'SC - HRF in V1 CF'; 'SC - HRF'; ...
    'SC - HRF difference'; 'SC - CF temporal offset'};
for hh = 1:length(hemis)
    ct = 0;
    for i = 1:4
        ct = ct + 1;
        data(:,ct) = datamat(:,i,hh,1);
        ct = ct + 1;
        data(:,ct) = datamat(:,i,hh,2);
        ct = ct + 1;
        data(:,ct) = datamat(:,i,hh,1) - datamat(:,i,hh,2);
        ct = ct + 1;
        data(:,ct) = datamat(:,i,hh,3);
    end
    f = figure;
    t = uitable(f, 'Data', data, 'ColumnName',xlabels,'RowName',rownames,'Position',[100 500 1720 100]);
end
%% Individual scatter plots
subjects = {'AEK' 'ASB' 'GKA'};
session_dirs = {...
    '/jet/abock/data/Retinotopy/AEK/10012014' ...
    '/jet/abock/data/Retinotopy/ASB/10272014' ...
    '/jet/abock/data/Retinotopy/GKA/06052015'};
ROIs = {'V2' 'V3' 'LGN' 'SC'};
lims = [-10 10];
for ss=1:length(subjects)
    figure;
    session_dir = session_dirs{ss};
    for rr=1:length(ROIs);
        ROI = ROIs{rr};
        if strcmp(ROI,'V2') || strcmp(ROI,'V3')
            map_type = 'cortex';
        elseif strcmp(ROI,'LGN') || strcmp(ROI,'SC')
            map_type = 'volume';
        end
        tmp = get_HRF_vals(session_dir,'lh',map_type,'prf_V1','movie',ROI);
        rh = get_HRF_vals(session_dir,'rh',map_type,'prf_V1','movie',ROI);
        X = [tmp.cfmask-tmp.prfmask;rh.cfmask - rh.prfmask];
        Y = [tmp.cfmask_shift;rh.cfmask_shift];
        P = polyfit(X,Y,1);
        pX = lims;
        pY = polyval(P,pX);
        subplot(2,2,rr);
        plot(X,Y,'.');hold on
        plot(pX,pY,'r');
        xlim(lims);
        ylim(lims);
        axis square;
        title(ROI);
        legend(['corr = ' num2str(corr(X,Y))],['slope = ' num2str(P(1))]);
    end
end
%% Subject average scatter plots
clear X Y pX pY
subjects = {'AEK' 'ASB' 'GKA'};
session_dirs = {...
    '/jet/abock/data/Retinotopy/AEK/10012014' ...
    '/jet/abock/data/Retinotopy/ASB/10272014' ...
    '/jet/abock/data/Retinotopy/GKA/06052015'};
ROIs = {'V2' 'V3' 'LGN' 'SC'};
hemis = {'lh' 'rh'};
for ss=1:length(subjects)
    session_dir = session_dirs{ss};
    ct = 0;
    for rr=1:length(ROIs);
        ROI = ROIs{rr};
        if strcmp(ROI,'V2') || strcmp(ROI,'V3')
            map_type = 'cortex';
        elseif strcmp(ROI,'LGN') || strcmp(ROI,'SC')
            map_type = 'volume';
        end
        for hh = 1:length(hemis)
            tmp = get_HRF_vals(session_dir,hemis{hh},map_type,'prf_V1','movie',ROI);
            dataX(ss,rr,hh) = median(tmp.cfmask-tmp.prfmask);
            dataY(ss,rr,hh) = median(tmp.cfmask_shift);
        end
    end
end
%%
lims=[-1 2];
mSize = 15;
markColors = {'b' 'r' 'k'};
markSyms = {'^','s','o','p'};
figure;
for ss = 1:length(subjects)
    for rr=1:length(ROIs); 
    plot(dataX(ss,rr,1),dataY(ss,rr,1),[markSyms{rr} markColors{ss}],'MarkerFaceColor',markColors{ss},'MarkerSize',mSize);hold on
    plot(dataX(ss,rr,2),dataY(ss,rr,2),[markSyms{rr} markColors{ss}],'MarkerSize',mSize);hold on
    end
end
xlim(lims);
ylim(lims);
axis square
legend(...
    'AEK - lh V2','AEK - rh V2','AEK - lh V3','AEK - rh V3',...
    'AEK - lh LGN','AEK - rh LGN','AEK - lh SC','AEK - rh SC', ...
    'ASB - lh V2','ASB - rh V2','ASB - lh V3','ASB - rh V3',...
    'ASB - lh LGN','ASB - rh LGN','ASB - lh SC','ASB - rh SC', ...
    'GKA - lh V2','GKA - rh V2','GKA - lh V3','GKA - rh V3',...
    'GKA - lh LGN','GKA - rh LGN','GKA - lh SC','GKA - rh SC'...
    );
P = polyfit(dataX(:),dataY(:),1);
pX = lims;
pY = polyval(P,pX);
plot(pX,pY,'r');
title(['Slope = ' num2str(P(1))]);
%%
pX = lims;
pY = polyval(P,pX);
plot(X,Y,'.');hold on
plot(pX,pY,'r');
xlim(lims);set(gca,'XTick',lims(1):1:lims(2));
ylim(lims);set(gca,'YTick',lims(1):1:lims(2));
axis square;
title(subjects{ss},'FontSize',20);
legend(['corr = ' num2str(corr(X,Y))],['slope = ' num2str(P(1))]);
xlabel('Mean peak HRF difference (V1 - ROI) (seconds)','FontSize',15);
ylabel('Mean CF temporal offset (seconds)','FontSize',15);
Xsum(ss,:) = X;
Ysum(ss,:) = Y;
columnNames = {'lh V2' 'rh V2' 'lh V3' 'rh V3' 'lh LGN' 'rh LGN' 'lh SC' 'rh SC'};
f = figure;
t = uitable(f, 'Data', Xsum, 'ColumnName',columnNames,'RowName',subjects,'Position',[585 500 750 100]);
f = figure;
t = uitable(f, 'Data', Ysum, 'ColumnName',columnNames,'RowName',subjects,'Position',[585 500 750 100]);


