function MultiplotWithError (output_dir,dataFiles, errorFiles, legendTexts, titleText, saveName)

%% set axes
figure('units','normalized','position',[0 0 1 1]);
hold on
xlabel('Time in seconds');
ylabel('Percent Signal Change');
title (titleText,'Interpreter','none')
xlims = [-1 15];
ylims = [-0.5 1.2];
xTick = [0 1 2 3 4 5 6 7 8 9 10 11 12 13];
xLabels = xTick;
h = plot([xlims(1) xlims(end)],[0 0],'k');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xlim(xlims); ylim(ylims);
ax = gca;
set(ax,'XTick',xTick);
set(ax,'XTickLabel',xLabels);

offsets = 0:(0.5/(length(dataFiles)-1)):0.5;

for ff = 1:length(dataFiles)
            y(:,ff) = csvread(dataFiles{ff});
            dataP = y(:,ff);
            e(:,ff) = csvread(errorFiles{ff});
            errorP = e(:,ff);
    end
    dataL = length(dataP);
    x = 0:1:dataL-1;
    g = errorbar(x+offsets(ff),dataP,errorP,'.k','MarkerSize',16);% hold on;
set(get(get(g,'Annotation'),'LegendInformation'),'IconDisplayStyle','children');
p = plot(x+offsets(ff), dataP, 'MarkerSize',16);
    axis square;
    legendInfo{ff} = (legendTexts{ff});

legend (p,legendInfo, 'Interpreter','none');

%save pdf in output dir

if ~exist (output_dir,'dir')
    mkdir (output_dir);
end

set(gcf, 'PaperPosition', [0 0 7 7]);
set(gcf, 'PaperSize', [7 7]);
saveas(gcf, fullfile(output_dir, saveName), 'pdf'); 
close all;