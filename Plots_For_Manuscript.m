%% Plot figures from the report:

%% Figure 1 - Single Summary statistic posteriors and bar plot with varying number of rows

PlotSinglePosterior(24,1,1);
PlotSinglePosterior(24,1,7);
PlotSinglePosterior(24,1,17);
PlotSinglePosterior(6,1,17);
PlotSinglePosterior(24,1,15);
PlotSinglePosterior(6,1,15);
PlotBarsVarySD(1);

%% Figure 2 - Vary Density

PlotSinglePosterior(24,1,17);
PlotSinglePosterior(24,3,17);
PlotBarsVaryDF(24);
PlotSinglePosterior(6,1,17);
PlotSinglePosterior(6,3,17);
PlotBarsVaryDF(6);

%% Figure 3 - Multiple SS posteriors

PlotSinglePosterior(24,1,[1,14]);
PlotSinglePosterior(24,1,[1,7]);
PlotSinglePosterior(24,1,[2,17]);
PlotSinglePosterior(24,1,[2,7]);
PlotBarsMaxIKL;

%% Figure 4 - ABC-DC

PlotABCDCchains;
PlotComparisonPosteriors;

%% Export code:
% For all of the above except figure 4 posteriors, export using:
addpath export_fig
export_fig( gcf, 'filename', '-painters','-nocrop','-transparent','-pdf' );

% For figure 4 posteriors, export using:
h=figure(NUM); % where NUM is the number of the figure you're saving
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'filename','-dpdf','-r0')