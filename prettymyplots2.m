function prettymyplots2 (Nx,Ny,prefX,prefY)

%Ny=7;
%Nx=7;
%prefX='% 1.1f';
%prefY='% 1.1f';

box on
set(gca,'TickLabelInterpreter','latex') %% NEW
set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'Color','none')
set(gca,'FontSize',20,'fontweight','normal')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','normal')
set(gca,'YTickMode','manual')
set(gca,'XTickMode','manual')
set(gcf,'color',[1 1 1])
set(gcf,'position',[1 40 560 420])
%%


yl = ylim;  
xl = xlim;

e=log10(yl(2));
e=floor(e);
if abs(e)<2
    e=0;
end
yt=get(gca,'ytick')/10^e;
set(gca,'ytick',linspace(yl(1),yl(2),Ny))
newyt=linspace(yt(1),yt(end),Ny);
ytl=cell([1,Ny]);
for j=1:length(ytl)
    ytl{j}=sprintf(prefY,newyt(j));
end
set(gca,'yticklabel',ytl);
fs = get(gca,'fontsize');
set(gca,'units','normalized');


xsl=diff(xl)/6;
if e~=0
    text(xl(1)-.1*xsl,yl(2),sprintf('$\\times10^{%d}$',e),'fontsize',fs,'VerticalAlignment','bottom','interpreter','latex');
end

%%

e=log10(xl(2));
e=floor(e);
if abs(e)<2
    e=0;
end
xt=get(gca,'xtick')/10^e;
set(gca,'xtick',linspace(xl(1),xl(2),Nx))
newxt=linspace(xt(1),xt(end),Nx);
xtl=cell([1,Nx]);
for j=1:length(xtl)
    xtl{j}=sprintf(prefX,newxt(j));
end
set(gca,'xticklabel',xtl);
fs = get(gca,'fontsize');
set(gca,'units','normalized');
yl = ylim;


ysl=diff(yl)/6;
if e~=0
    text(xl(2)-xsl/2,yl(1)-1.4*ysl,sprintf('$\\times10^{%d}$',e),'fontsize',fs,'VerticalAlignment','bottom','interpreter','latex');
end