function MultiplotWithError(output_dir,dataFiles, errorFiles, legendTexts, titleText, saveName)
% MultiplotWithError(output_dir,dataFiles, errorFiles, legendTexts, titleText, saveName)
%
% <GF> Please add input arguments and usage here
%
% Input arguments:
% ================
%
%   session_dir : 
%   ...
%
% Usage:
% ======
%
% <GF>
%
%
% 7/2/2016  gf, ms      Written and commented.

%% set axes
firFig = figure('units','normalized','position',[0 0 1 1]);
hold on;
xlabel('Time [sec]');
ylabel('Signal change [%]');
title(titleText,'Interpreter','none')
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
legendInfo{ff} = (legendTexts{ff});
legend (p,legendInfo, 'Interpreter','none');

% Save pdf in output dir
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% Adjust the figure
adjustPlot(firFig);
saveas(firFig, fullfile(output_dir, saveName), 'pdf');
close all;