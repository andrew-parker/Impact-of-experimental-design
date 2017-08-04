function prettymyplots
%title(' ');
set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'Color','none')
fs=20;
set(gca,'FontSize',fs,'fontweight','normal')
set(findall(gcf,'type','text'),'FontSize',fs,'fontWeight','normal')
set(gca,'YTickMode','manual')
set(gca,'XTickMode','manual')
%% set tick labels
yl=ylim;  
% get order of magnitude
e=log10(yl(2));
%e=sign(e)*floor(abs(e));
e=floor(e);
% get and rescale yticks
yt=get(gca,'ytick')/10^e;
% create tick labels
Ny=6;
set(gca,'ytick',linspace(yl(1),yl(2),Ny))
newyt=linspace(yt(1),yt(end),Ny);
ytl=cell([1,Ny]);
for j=1:length(ytl)
% the space after the percent gives the same size to positive and
% negative numbers. The number of decimal digits can be changed.
ytl{j}=sprintf('% 1.1f',newyt(j));
end
% set tick labels
set(gca,'yticklabel',ytl);
% place order of magnitude
set(gca,'fontsize',fs);
set(gca,'units','normalized');
xl = xlim;
%text(xl(1),yl(2)+0.0001,sprintf('$\\times10^{%d}$',e),'fontsize',fs,'VerticalAlignment','bottom','interpreter','latex');
text(xl(1)-0.03,yl(2)+0.0002,sprintf('$\\times10^{%d}$',e),'fontsize',fs,'VerticalAlignment','bottom','interpreter','latex');
%% set xtick labels 0.0, 0.5, 1.0
Nx=3;
xt=get(gca,'xtick');
if (xl(2)==0.1)||(xl(2)==0.5)
set(gca,'xtick',linspace(xt(1),xt(end),Nx))
set(gca,'xticklabel',num2str(get(gca,'xtick')','%1.2f'))
elseif (xl(2)==0.05)

e=log10(xl(2));
%e=sign(e)*floor(abs(e));
e=floor(e);
% get and rescale yticks
xt=get(gca,'xtick')/10^e;
% create tick labels
set(gca,'xtick',linspace(xl(1),xl(2),Nx))
newxt=linspace(xt(1),xt(end),Nx);
xtl=cell([1,Nx]);
for j=1:length(xtl)
% the space after the percent gives the same size to positive and
% negative numbers. The number of decimal digits can be changed.
xtl{j}=sprintf('% 1.1f',newxt(j));
end
% set tick labels
set(gca,'xticklabel',xtl);
% place order of magnitude
set(gca,'fontsize',fs);
set(gca,'units','normalized');

text(xl(2),-0.2*yl(2),sprintf('$\\times10^{%d}$',e),'fontsize',fs,'VerticalAlignment','bottom','interpreter','latex');
%text(xl(2)*1,-yl(2)*0.15,sprintf('\\times10^{-2}'),'FontSize',20)

else
set(gca,'xtick',linspace(xt(1),xt(end),Nx))
set(gca,'xticklabel',num2str(get(gca,'xtick')','%1.1f'))
end
%% set colorbar
hcb=colorbar;
yla=get(hcb,'YTick');
ylb=[yla(1),yla(end)];
% get order of magnitude
e=log10(ylb(2));
e=sign(e)*floor(abs(e));
% get and rescale yticks
yt=get(hcb,'ytick')/10^e;
% create tick labels
ytl=cell(size(yt));
for j=1:length(yt)
% the space after the percent gives the same size to positive and
% negative numbers. The number of decimal digits can be changed.
ytl{j}=sprintf('% 1.1f',yt(j));
end
% set tick labels
set(hcb,'yticklabel',ytl);
set(hcb,'ticklabelinterpreter','latex')
% place order of magnitude
set(hcb,'fontsize',fs);
set(hcb,'units','normalized');
%xl = xlim;
%text(xl(1),yl(2),sprintf('\\times10^{%d}',e),'fontsize',fs,'VerticalAlignment','bottom');
%text(xl(2)*1.15,yl(2)*1.01,sprintf('$\\times10^{%d}$',e),'fontsize',fs,'VerticalAlignment','bottom','interpreter','latex');

text(xl(2)*1.03,yl(2)*1.01,sprintf('$\\times10^{%d}$',e),'fontsize',fs,'VerticalAlignment','bottom','interpreter','latex');

set(gca,'FontSize',fs,'fontweight','normal')
set(findall(gcf,'type','text'),'FontSize',fs,'fontWeight','normal')