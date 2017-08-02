%% Plot accepted parameter chains from ABC-DC

load('ABCsimulationsCXYSD6DF1.mat', 'thetaRecord')

len=10000;
figure()
plot(thetaRecord(1,:)')
set(gcf,'position',[1 40 680 420]) %changed to 680 to fit 6 ticks
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','normal')
set(gca,'FontSize',20,'fontWeight','normal')
set(gca,'TickLabelInterpreter', 'latex');
ylabel('$P_m$','interpreter','latex')
%xlabel('Iteration','interpreter','latex')
hold on
plot(1:len,0.25*ones(1,len),'--','linewidth',2)%,'color',0.8*[1 1 1])
plot(5000*ones(1,1001),0:0.001:1,'-','color',0*[1 1 1],'linewidth',2)
%plot(14000*ones(1,1001),0:0.001:1,'-','color',0*[1 1 1])
set(gca,'yticklabels',num2str(get(gca,'ytick')','%1.1f'));
%set(gca,'xticklabels',{'0.0','0.5','1.0','1.5','2.0','2.5'})
set(gca,'xticklabels',{'0','2','4','6','8','10'})
text(len,-.0036, '$\times 10^{3}$','interpreter','latex')
%text(2500,1, '$\arrowdown \delta$','interpreter','latex')
text(2500,1.04, '$\downarrow \delta$','interpreter','latex','fontsize',20)
text(7500,1.04, '$\uparrow K$','interpreter','latex','fontsize',20)
set(gca,'yticklabels',num2str(get(gca,'ytick')','%1.1f'));
%grid on
set(gcf,'color',[1 1 1])
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','normal')


figure()
plot(100*thetaRecord(2,:)')
set(gcf,'position',[1 40 680 420]) %changed to 680 to fit 6 ticks
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','normal')
set(gca,'FontSize',20,'fontWeight','normal')
set(gca,'TickLabelInterpreter', 'latex');
ylabel('$P_p$','interpreter','latex')
%xlabel('Iteration','interpreter','latex')
hold on
plot(1:len,0.25*ones(1,len),'--','linewidth',2)%,'color',0.8*[1 1 1])
plot(5000*ones(1,1001),0:0.001:1,'-','color',0*[1 1 1],'linewidth',2)

set(gca,'yticklabels',num2str(get(gca,'ytick')','%1.1f'));

set(gca,'xticklabels',{'0','2','4','6','8','10'})
text(len,-.0036, '$\times 10^{3}$','interpreter','latex')

set(gca,'yticklabels',num2str(get(gca,'ytick')','%1.1f'));
%grid on
set(gcf,'color',[1 1 1])
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','normal')
text(0,1.033, '$\times 10^{-2}$','interpreter','latex','fontsize',20)
