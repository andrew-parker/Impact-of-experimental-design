function PlotBarsVarySD(DF)


% information1RE is copied into info6 and info 24
% info6=information1RE; %e.t.c.
load(['infoRE_SD24DF',int2str(DF),'.mat'])
info24=information1RE;
load(['infoRE_SD12DF',int2str(DF),'.mat'])
info12=information1RE;
load(['infoRE_SD6DF',int2str(DF),'.mat'])
info6=information1RE;
%infomat=[mean(info24([1 7 15 17],1,:),3),mean(info6([1 7 15 17],1,:),3)];
%infomatstd=[std(info24([1 7 15 17],1,:),[],3),std(info6([1 7 15 17],1,:),[],3)];
symbols={'$N$';'$C_{XY}$';'\max{C}';'$\hat{C}$';'$\lambda$';'\lambda_2';'$\|x\|$';'\|x\|_2$';'$\Gamma$';'$\Gamma_2$';'$\kappa_4$';'$\kappa_8$';'$Q_2$';'$Q_4$';'$Q_8$';'$\hat{C_{XY}}$';'$C_Y$';'$C_X$';'$\hat{C_Y}$';'$\hat{C_X}$'};
vectorSS=[1 7 15 17];
vectorSS=[1,11:15,7,9,5,2,16,17,19]'; %no 3,4,6,8,10,18,20

infomat=[mean(info24(vectorSS,1,:),3),mean(info12(vectorSS,1,:),3),mean(info6(vectorSS,1,:),3)];
infomatstd=[std(info24(vectorSS,1,:),[],3),std(info12(vectorSS,1,:),[],3),std(info6(vectorSS,1,:),[],3)];

figure()
hbar=bar(infomat);
ylim([0 4])
set(gca,'xticklabel',symbols(vectorSS))
set(gca,'TickLabelInterpreter','latex') %% NEW
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.1f'))
grid on
colormap summer
cmap=colormap;
set(hbar,{'FaceColor'},{cmap(1,:);cmap(32,:);cmap(64,:)});
hold on
%set(gcf,'position',[1 40 560 420])
% errorbar([-1,1:length(vectorSS),6]-.22,[-3;infomat(:,1);-3],[0;infomatstd(:,1);0],'linestyle','none','color',.1*[1 1 1],'linewidth',2)
% errorbar([-1,1:length(vectorSS),6],[-3;infomat(:,2);-3],[0;infomatstd(:,2);0],'linestyle','none','color',.1*[1 1 1],'linewidth',2)
% errorbar([-1,1:length(vectorSS),6]+.22,[-3;infomat(:,3);-3],[0;infomatstd(:,3);0],'linestyle','none','color',.1*[1 1 1],'linewidth',2)
set(gcf,'position',[1 40 1120 420])
errorbar([1:length(vectorSS)]-.22,[infomat(:,1)],[infomatstd(:,1)],'linestyle','none','color',.1*[1 1 1],'linewidth',2)
errorbar([1:length(vectorSS)],[infomat(:,2)],[infomatstd(:,2)],'linestyle','none','color',.1*[1 1 1],'linewidth',2)
errorbar([1:length(vectorSS)]+.22,[infomat(:,3)],[infomatstd(:,3)],'linestyle','none','color',.1*[1 1 1],'linewidth',2)
%errorbar([-2,1:4,7]+.15,[-3;infomat(:,2);-3],[0;infomatstd(:,2);0],'linestyle','none','color',cmap(1,:),'linewidth',2)
xlim([0.5, length(vectorSS)+.5]);
set(gcf,'color',[1 1 1])
ylabel('$I_{KL}$','interpreter','latex')


%% Editting
%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'Color','none')
fs=20;
set(gca,'FontSize',fs,'fontweight','normal')
set(findall(gcf,'type','text'),'FontSize',fs,'fontWeight','normal')
set(gca,'YTickMode','manual')
set(gca,'XTickMode','manual')

yl=ylim;  
% get order of magnitude
e=log10(yl(2));
%e=sign(e)*floor(abs(e));
e=floor(e);
% get and rescale yticks
yt=get(gca,'ytick')/10^e;
% create tick labels
Ny=5;
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
