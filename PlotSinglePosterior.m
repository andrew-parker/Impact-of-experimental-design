function PlotSinglePosterior(SD,DF,Indices)

symbols={'$N$';'$C_{XY}$';'\max{C}';'$\hat{C}$';'$\lambda$';'\lambda_2';'$\|x\|$';'\|x\|_2$';'$\tau$';'$\tau_2$';'$\kappa_4$';'$\kappa_8$';'$Q_2$';'$Q_4$';'$Q_8$';'$\hat{C_{XY}}$';'$C_Y$';'$C_X$';'$\hat{C_Y}$';'$\hat{C_X}$'};
%COMBOS1=nchoosek(1:20,1);
COMBOS2=nchoosek(1:20,2);
COMBOS3=nchoosek(1:20,3);

numstats=length(Indices);
if numstats==1
    Index=Indices;
elseif numstats==2
    Index=intersect(find(COMBOS2(:,1)==Indices(1)),find(COMBOS2(:,2)==Indices(2)));
elseif numstats==3
    Index=intersect(find(COMBOS3(:,1)==Indices(1)),find(COMBOS3(:,2)==Indices(2)),find(COMBOS3(:,3)==Indices(3)));
end

load(['infoRE_SD',num2str(SD),'DF',num2str(DF),'.mat']);

%% Figure 1 
%numstats=1;

%load('infoRE_SD24DF1.mat')
%load('infoRE_SD6DF1.mat')
%Index=1;
%Index=7;
%Index=17;
%Index=15;

%% Figure 2
%numstats=1;
%Index=17;
%load('infoRE_SD24DF3.mat')
%load('infoRE_SD6DF3.mat')

%% Figure 3
%numstats=2;
%load('infoRE_SD24DF1.mat')
%Index=intersect(find(COMBOS2(:,1)==1),find(COMBOS2(:,2)==14));
%Index=intersect(find(COMBOS2(:,1)==1),find(COMBOS2(:,2)==7));
%load('infoRE_SD6DF1.mat')
%Index=intersect(find(COMBOS2(:,1)==2),find(COMBOS2(:,2)==17));
%Index=intersect(find(COMBOS2(:,1)==2),find(COMBOS2(:,2)==7));



%%

Dlarge=mean(DlargeRE,3);


%accCI1a=[1,11:15,7,9,5,2,16,17,19]'; %no 3,4,6,8,10,18,20
%accCI2a=find(ones(length(COMBOS2),1)-(COMBOS2(:,1)==3 | COMBOS2(:,1)==4 | COMBOS2(:,1)==6 | COMBOS2(:,1)==8 | COMBOS2(:,1)==10 | COMBOS2(:,1)==18 | COMBOS2(:,1)==20 | COMBOS2(:,2)==3 | COMBOS2(:,2)==4 | COMBOS2(:,2)==6 | COMBOS2(:,2)==8 | COMBOS2(:,2)==10 | COMBOS2(:,2)==18 | COMBOS2(:,2)==20));
%accCI3a=find(ones(length(COMBOS3),1)-(COMBOS3(:,1)==3 | COMBOS3(:,1)==4 | COMBOS3(:,1)==6 | COMBOS3(:,1)==8 | COMBOS3(:,1)==10 | COMBOS3(:,1)==18 | COMBOS3(:,1)==20 | COMBOS3(:,2)==3 | COMBOS3(:,2)==4 | COMBOS3(:,2)==6 | COMBOS3(:,2)==8 | COMBOS3(:,2)==10 | COMBOS3(:,2)==18 | COMBOS3(:,2)==20 | COMBOS3(:,3)==3 | COMBOS3(:,3)==4 | COMBOS3(:,3)==6 | COMBOS3(:,3)==8 | COMBOS3(:,3)==10 | COMBOS3(:,3)==18 | COMBOS3(:,3)==20 ));


if numstats==1
    STC=Index;
elseif numstats==2
    STC=COMBOS2(Index,:);
elseif numstats==3
    STC=COMBOS3(Index,:);
end

%load('infoRE_SD6DF1.mat')




Ratio=100;
P1=0.25;
P2=0.0025;
P1min=0;
P1max=1;%0.5
P2min=0;
P2max=0.01;%0.005
xlab='$P_m$';
ylab='$P_p$';
method=1;


%figure(2)
set(gcf,'position',[1 40 560 420])
%title('Posteriors computed with one summary statistic')
%fignum=0;
%STC=17;
%for i=STC
  %  figure()
   % figure(fignum)
%    fignum=fignum+1;
    %[~,I]=sort(DlargeRE(j,:,i));
%     switch numstats
%         case 1
%             STC=COMBOS1(Index,:);
%             Stitle=strjoin(symbols(STC)');
%         case 2
%             STC=COMBOS2(Index,:);
%             Stitle=strjoin(symbols(COMBOS2(Index,:))');
%         case 3
%             STC=COMBOS3(Index,:);
%             Stitle=strjoin(symbols(COMBOS3(Index,:))');
%     end
    if numstats==1
        ca2=symbols(STC)';
    end
    if numstats==2
        ca=symbols(STC)';
        ca2={ca{1},'\&',ca{2}};
    end
    if numstats==3
        ca=symbols(STC)';
        ca2={ca{1},'\&',ca{2},'\&',ca{3}};
    end
    %Stitle=strjoin(symbols(STC)');
    Stitle=strjoin(ca2);
    D=sum(Dlarge(STC,:).^2,1);
    [~,I]=sort(D);
        
    k=200;
    data=[theta(1,I(1:k))',theta(2,I(1:k))'];
    if method==1
        [bandwidth,density,X,Y]=kde2d(data,512,[P1min,P2min],[P1max,P2max]);
        contourf(X,Y,density,'LineStyle','None');
    else
        %data=[theta(1,I(1:k))',theta(2,I(1:k))'];
        p=20; %the resolution of the discretisation
        density=zeros(p);
        for j=1:k
           % density(ceil(p*data(j,1)),ceil(p*100*data(j,2)))=density(ceil(p*data(j,1)),ceil(p*100*data(j,2)))+1; %Pm,Pp
            %density(ceil(p*data(j,1)),ceil(p*(data(j,2)+0.75)))=density(ceil(p*data(j,1)),ceil(p*(data(j,2)+0.75)))+1;
            density(ceil(p*(data(j,1)-P1min)/(P1max-P1min)),ceil(p*(data(j,2)-P2min)/(P2max-P2min)))=density(ceil(p*(data(j,1)-P1min)/(P1max-P1min)),ceil(p*(data(j,2)-P2min)/(P2max-P2min)))+1;
        end
        contourf(linspace(P1min,P1max,p),linspace(P2min,P2max,p),density','LineStyle','None');
    end
    colormap summer, hold on, alpha(.8)
    cmap=colormap;
    %set(gca, 'color', 'blue');
    set(gca,'color',cmap(1,:))
    plot(data(:,1),data(:,2),'w.','MarkerSize',3)
    plot(0.25*ones(512,1),linspace(0,0.01,512),'color',[.8 .8 .8],'linestyle','--','linewidth',1.5)
    plot(linspace(0,1,512),0.0025*ones(512,1),'color',.8*[1 1 1],'linestyle','--','linewidth',1.5)
    %plot(P1,P2,'r*','MarkerSize',5)
    xlim([P1min,P1max])
    ylim([P2min,P2max])
    %ylim([0,0.005])
    xlabel(xlab,'interpreter','latex')
    ylabel(ylab,'interpreter','latex')
    title(Stitle,'interpreter','latex')
    set(gca,'TickLabelInterpreter','latex') %% NEW
    colorbar
    grid on
    pause(0.1);
    hold off    
    prettymyplots; %% NEW
%end
set(gcf,'color',[1 1 1])
set(gcf,'position',[1 40 560 420])