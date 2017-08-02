
%h=bar([1 2 3; 1 2 3;1 2 3]);
%cmapshort=get(gca,'colororder')

h=bar([1 2 3; 1 2 3;1 2 3]),
colormap(cmapshort(1:3,1:3));

%colormap summer
%colormap cmapshort
%cmap=colormap;
%cmap2=0.5*[1 1 1]+0.5*cmap([1;32;64],:);
%colormap(cmap2)
%hl=legend(' 24 rows',' 12 rows',' 6 rows');
%hl=legend('$\,$24 rows','$\,$12 rows','$\,$6 rows');
hl=legend('$\,$ABC-rej','$\,$ABC-DC, K=1','$\,$ABC-DC, K=4');
%hl=legend('$\,$24 cells','$\,$48 cells','$\,$72 cells');
set(gca,'position',[0 0 1 0.0])
set(hl,'interpreter','latex','box','off','position',[0.0 0.4 1 .3],'fontsize',24);
set(gcf,'color',[1 1 1])

%addpath export_fig
%export_fig( gcf, ['legendattempt'], '-painters','-nocrop','-transparent','-pdf' );