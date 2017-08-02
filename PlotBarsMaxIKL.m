%% Plot best combination IKLs for each experiment:
%for Experiment=1:3
%Experiment=3;
Experiment=1;
if Experiment==1
AllInformation1=zeros(20,3,3);
AllInformation2=zeros(190,3,3);
AllInformation3=zeros(1140,3,3);
AllInformation1std=zeros(20,3,3);
AllInformation2std=zeros(190,3,3);
AllInformation3std=zeros(1140,3,3);
load('infoRE_SD24DF1.mat')
AllInformation1(:,:,1)=mean(information1RE,3);
AllInformation2(:,:,1)=mean(information2RE,3);
AllInformation3(:,:,1)=mean(information3RE,3);
AllInformation1std(:,:,1)=std(information1RE,[],3);
AllInformation2std(:,:,1)=std(information2RE,[],3);
AllInformation3std(:,:,1)=std(information3RE,[],3);
load('infoRE_SD12DF1.mat')
AllInformation1(:,:,2)=mean(information1RE,3);
AllInformation2(:,:,2)=mean(information2RE,3);
AllInformation3(:,:,2)=mean(information3RE,3);
AllInformation1std(:,:,2)=std(information1RE,[],3);
AllInformation2std(:,:,2)=std(information2RE,[],3);
AllInformation3std(:,:,2)=std(information3RE,[],3);
load('infoRE_SD6DF1.mat')
AllInformation1(:,:,3)=mean(information1RE,3);
AllInformation2(:,:,3)=mean(information2RE,3);
AllInformation3(:,:,3)=mean(information3RE,3);
AllInformation1std(:,:,3)=std(information1RE,[],3);
AllInformation2std(:,:,3)=std(information2RE,[],3);
AllInformation3std(:,:,3)=std(information3RE,[],3);
% load('infoRE_SD3DF1.mat')
% AllInformation1(:,:,4)=mean(information1RE,3);
% AllInformation2(:,:,4)=mean(information2RE,3);
% AllInformation3(:,:,4)=mean(information3RE,3);
% AllInformation1std(:,:,4)=std(information1RE,[],3);
% AllInformation2std(:,:,4)=std(information2RE,[],3);
% AllInformation3std(:,:,4)=std(information3RE,[],3);
elseif Experiment==2
AllInformation1=zeros(20,3,3);
AllInformation2=zeros(190,3,3);
AllInformation3=zeros(1140,3,3);
AllInformation1std=zeros(20,3,3);
AllInformation2std=zeros(190,3,3);
AllInformation3std=zeros(1140,3,3);
load('infoRE_SD24DF1.mat')
AllInformation1(:,:,1)=mean(information1RE,3);
AllInformation2(:,:,1)=mean(information2RE,3);
AllInformation3(:,:,1)=mean(information3RE,3);
AllInformation1std(:,:,1)=std(information1RE,[],3);
AllInformation2std(:,:,1)=std(information2RE,[],3);
AllInformation3std(:,:,1)=std(information3RE,[],3);
load('infoRE_SD24DF2.mat')
AllInformation1(:,:,2)=mean(information1RE,3);
AllInformation2(:,:,2)=mean(information2RE,3);
AllInformation3(:,:,2)=mean(information3RE,3);
AllInformation1std(:,:,2)=std(information1RE,[],3);
AllInformation2std(:,:,2)=std(information2RE,[],3);
AllInformation3std(:,:,2)=std(information3RE,[],3);
load('infoRE_SD24DF3.mat')
AllInformation1(:,:,3)=mean(information1RE,3);
AllInformation2(:,:,3)=mean(information2RE,3);
AllInformation3(:,:,3)=mean(information3RE,3);
AllInformation1std(:,:,3)=std(information1RE,[],3);
AllInformation2std(:,:,3)=std(information2RE,[],3);
AllInformation3std(:,:,3)=std(information3RE,[],3);
else
AllInformation1=zeros(20,3,3);
AllInformation2=zeros(190,3,3);
AllInformation3=zeros(1140,3,3);
AllInformation1std=zeros(20,3,3);
AllInformation2std=zeros(190,3,3);
AllInformation3std=zeros(1140,3,3);
load('infoRE_SD6DF1.mat')
AllInformation1(:,:,1)=mean(information1RE,3);
AllInformation2(:,:,1)=mean(information2RE,3);
AllInformation3(:,:,1)=mean(information3RE,3);
AllInformation1std(:,:,1)=std(information1RE,[],3);
AllInformation2std(:,:,1)=std(information2RE,[],3);
AllInformation3std(:,:,1)=std(information3RE,[],3);
load('infoRE_SD6DF2.mat')
AllInformation1(:,:,2)=mean(information1RE,3);
AllInformation2(:,:,2)=mean(information2RE,3);
AllInformation3(:,:,2)=mean(information3RE,3);
AllInformation1std(:,:,2)=std(information1RE,[],3);
AllInformation2std(:,:,2)=std(information2RE,[],3);
AllInformation3std(:,:,2)=std(information3RE,[],3);
load('infoRE_SD6DF3.mat')
AllInformation1(:,:,3)=mean(information1RE,3);
AllInformation2(:,:,3)=mean(information2RE,3);
AllInformation3(:,:,3)=mean(information3RE,3);
AllInformation1std(:,:,3)=std(information1RE,[],3);
AllInformation2std(:,:,3)=std(information2RE,[],3);
AllInformation3std(:,:,3)=std(information3RE,[],3);
end

%% Copy and paste the below, loading different information each time, to obtain the table in the report

%information1=AllInformation(newSTC,:,1);
%nformation2=
%accepted COMBOS2 Indices (not excluded)

COMBOS1=nchoosek(1:20,1);
COMBOS2=nchoosek(1:20,2);
COMBOS3=nchoosek(1:20,3);
symbols={'$N$';'$C_{XY}$';'\max{C}';'$\hat{C}$';'$\lambda$';'\lambda_2';'$\|x\|$';'\|x\|_2$';'$\Gamma$';'$\Gamma_2$';'$\kappa_4$';'$\kappa_8$';'$Q_2$';'$Q_4$';'$Q_8$';'$\hat{C_{XY}}$';'$C_Y$';'$C_X$';'$\hat{C_Y}$';'$\hat{C_X}$'};
AllMaxIK=[];
AllMaxIKa=[];
AllMaxSTDs=[];
AllMaxSTDsa=[];
IK2D=[];
IK2Da=[];
IKSTD2D=[];
IKSTD2Da=[];

for expSample=1:size(AllInformation1,3)
    
information1=AllInformation1(:,:,expSample);
information2=AllInformation2(:,:,expSample);
information3=AllInformation3(:,:,expSample);

%% Compute best combinations with all stats:
accCI1a=[1,11:15,7,9,5,2,16,17,19]'; %no 3,4,6,8,10,18,20
accCI2a=find(ones(length(COMBOS2),1)-(COMBOS2(:,1)==3 | COMBOS2(:,1)==4 | COMBOS2(:,1)==6 | COMBOS2(:,1)==8 | COMBOS2(:,1)==10 | COMBOS2(:,1)==18 | COMBOS2(:,1)==20 | COMBOS2(:,2)==3 | COMBOS2(:,2)==4 | COMBOS2(:,2)==6 | COMBOS2(:,2)==8 | COMBOS2(:,2)==10 | COMBOS2(:,2)==18 | COMBOS2(:,2)==20));
accCI3a=find(ones(length(COMBOS3),1)-(COMBOS3(:,1)==3 | COMBOS3(:,1)==4 | COMBOS3(:,1)==6 | COMBOS3(:,1)==8 | COMBOS3(:,1)==10 | COMBOS3(:,1)==18 | COMBOS3(:,1)==20 | COMBOS3(:,2)==3 | COMBOS3(:,2)==4 | COMBOS3(:,2)==6 | COMBOS3(:,2)==8 | COMBOS3(:,2)==10 | COMBOS3(:,2)==18 | COMBOS3(:,2)==20 | COMBOS3(:,3)==3 | COMBOS3(:,3)==4 | COMBOS3(:,3)==6 | COMBOS3(:,3)==8 | COMBOS3(:,3)==10 | COMBOS3(:,3)==18 | COMBOS3(:,3)==20 ));

[~,Ii1aa]=sort(information1(accCI1a,1));
[~,Ii1ba]=sort(information1(accCI1a,2));
[~,Ii1ca]=sort(information1(accCI1a,3));

[~,Ii2aa]=sort(information2(accCI2a,1));
[~,Ii2ba]=sort(information2(accCI2a,2));
[~,Ii2ca]=sort(information2(accCI2a,3));

[~,Ii3aa]=sort(information3(accCI3a,1));
[~,Ii3ba]=sort(information3(accCI3a,2));
[~,Ii3ca]=sort(information3(accCI3a,3));
TopK=0; %Top 3 <-> TopK=2, Top 1 <-> TopK=0;
TopSymbols1a=[symbols(accCI1a(Ii1aa(end-TopK:end)))';symbols(accCI1a(Ii1ba(end-TopK:end)))';symbols(accCI1a(Ii1ca(end-TopK:end)))'];
TopIK1a=[information1(accCI1a(Ii1aa(end-TopK:end)),1)';information1(accCI1a(Ii1ba(end-TopK:end)),2)';information1(accCI1a(Ii1ca(end-TopK:end)),3)'];
TopSTD1a=[AllInformation1std(accCI1a(Ii1aa(end-TopK:end)),1,expSample)';AllInformation1std(accCI1a(Ii1ba(end-TopK:end)),2,expSample)';AllInformation1std(accCI1a(Ii1ca(end-TopK:end)),3,expSample)'];

TopSymbols2a=[symbols(COMBOS2(accCI2a(Ii2aa(end-TopK:end)),:))';symbols(COMBOS2(accCI2a(Ii2ba(end-TopK:end)),:))';symbols(COMBOS2(accCI2a(Ii2ca(end-TopK:end)),:))'];
TopIK2a=[information2(accCI2a(Ii2aa(end-TopK:end)),1)';information2(accCI2a(Ii2ba(end-TopK:end)),2)';information2(accCI2a(Ii2ca(end-TopK:end)),3)'];
TopSTD2a=[AllInformation2std(accCI2a(Ii2aa(end-TopK:end)),1,expSample)';AllInformation2std(accCI2a(Ii2ba(end-TopK:end)),2,expSample)';AllInformation2std(accCI2a(Ii2ca(end-TopK:end)),3,expSample)'];

TopSymbols3a=[symbols(COMBOS3(accCI3a(Ii3aa(end-TopK:end)),:))';symbols(COMBOS3(accCI3a(Ii3ba(end-TopK:end)),:))';symbols(COMBOS3(accCI3a(Ii3ca(end-TopK:end)),:))'];
TopIK3a=[information3(accCI3a(Ii3aa(end-TopK:end)),1)';information3(accCI3a(Ii3ba(end-TopK:end)),2)';information3(accCI3a(Ii3ca(end-TopK:end)),3)'];
TopSTD3a=[AllInformation3std(accCI3a(Ii3aa(end-TopK:end)),1,expSample)';AllInformation3std(accCI3a(Ii3ba(end-TopK:end)),2,expSample)';AllInformation3std(accCI3a(Ii3ca(end-TopK:end)),3,expSample)'];
%% Compute best combinations excluding trajectory statistics:
accCI1=[1,11:15,5,2,16,17,19]'; %no 3,4,6,8,10,18,20 or 7 or 9
accCI2=find(ones(length(COMBOS2),1)-(COMBOS2(:,1)==3 | COMBOS2(:,1)==4 | COMBOS2(:,1)==6 | COMBOS2(:,1)==7 | COMBOS2(:,1)==8 | COMBOS2(:,1)==9 | COMBOS2(:,1)==10 | COMBOS2(:,1)==18 | COMBOS2(:,1)==20 | COMBOS2(:,2)==3 | COMBOS2(:,2)==4 | COMBOS2(:,2)==6 | COMBOS2(:,2)==7 | COMBOS2(:,2)==8 | COMBOS2(:,2)==9 | COMBOS2(:,2)==10 | COMBOS2(:,2)==18 | COMBOS2(:,2)==20));
accCI3=find(ones(length(COMBOS3),1)-(COMBOS3(:,1)==3 | COMBOS3(:,1)==4 | COMBOS3(:,1)==6 | COMBOS3(:,1)==7 | COMBOS3(:,1)==8 | COMBOS3(:,1)==9 | COMBOS3(:,1)==10 | COMBOS3(:,1)==18 | COMBOS3(:,1)==20 | COMBOS3(:,2)==3 | COMBOS3(:,2)==4 | COMBOS3(:,2)==6 | COMBOS3(:,2)==7 | COMBOS3(:,2)==8 | COMBOS3(:,2)==9 | COMBOS3(:,2)==10 | COMBOS3(:,2)==18 | COMBOS3(:,2)==20 | COMBOS3(:,3)==3 | COMBOS3(:,3)==4 | COMBOS3(:,3)==6 | COMBOS3(:,3)==7 | COMBOS3(:,3)==8 | COMBOS3(:,3)==9 | COMBOS3(:,3)==10 | COMBOS3(:,3)==18 | COMBOS3(:,3)==20 ));

[~,Ii1a]=sort(information1(accCI1,1));
[~,Ii1b]=sort(information1(accCI1,2));
[~,Ii1c]=sort(information1(accCI1,3));

[~,Ii2a]=sort(information2(accCI2,1));
[~,Ii2b]=sort(information2(accCI2,2));
[~,Ii2c]=sort(information2(accCI2,3));

[~,Ii3a]=sort(information3(accCI3,1));
[~,Ii3b]=sort(information3(accCI3,2));
[~,Ii3c]=sort(information3(accCI3,3));
TopK=0; %Top 3 <-> TopK=2, Top 1 <-> TopK=0;
TopSymbols1=[symbols(accCI1(Ii1a(end-TopK:end)))';symbols(accCI1(Ii1b(end-TopK:end)))';symbols(accCI1(Ii1c(end-TopK:end)))'];
TopIK1=[information1(accCI1(Ii1a(end-TopK:end)),1)';information1(accCI1(Ii1b(end-TopK:end)),2)';information1(accCI1(Ii1c(end-TopK:end)),3)'];
TopSTD1=[AllInformation1std(accCI1(Ii1a(end-TopK:end)),1,expSample)';AllInformation1std(accCI1(Ii1b(end-TopK:end)),2,expSample)';AllInformation1std(accCI1(Ii1c(end-TopK:end)),3,expSample)'];

TopSymbols2=[symbols(COMBOS2(accCI2(Ii2a(end-TopK:end)),:))';symbols(COMBOS2(accCI2(Ii2b(end-TopK:end)),:))';symbols(COMBOS2(accCI2(Ii2c(end-TopK:end)),:))'];
TopIK2=[information2(accCI2(Ii2a(end-TopK:end)),1)';information2(accCI2(Ii2b(end-TopK:end)),2)';information2(accCI2(Ii2c(end-TopK:end)),3)'];
TopSTD2=[AllInformation2std(accCI2(Ii2a(end-TopK:end)),1,expSample)';AllInformation2std(accCI2(Ii2b(end-TopK:end)),2,expSample)';AllInformation2std(accCI2(Ii2c(end-TopK:end)),3,expSample)'];

TopSymbols3=[symbols(COMBOS3(accCI3(Ii3a(end-TopK:end)),:))';symbols(COMBOS3(accCI3(Ii3b(end-TopK:end)),:))';symbols(COMBOS3(accCI3(Ii3c(end-TopK:end)),:))'];
TopIK3=[information3(accCI3(Ii3a(end-TopK:end)),1)';information3(accCI3(Ii3b(end-TopK:end)),2)';information3(accCI3(Ii3c(end-TopK:end)),3)'];
TopSTD3=[AllInformation3std(accCI3(Ii3a(end-TopK:end)),1,expSample)';AllInformation3std(accCI3(Ii3b(end-TopK:end)),2,expSample)';AllInformation3std(accCI3(Ii3c(end-TopK:end)),3,expSample)'];
%% Construct array from data:
% Should look like Best combination | IK | Best combo without traj ss | IK
% with three rows, for Pm,Pp, Pm and Pp
TopIKa=[TopIK1a,TopIK2a,TopIK3a]';
TopIndices(:,:,expSample)=[accCI1(Ii1a(end)),accCI1(Ii1b(end)),accCI1(Ii1c(end));accCI2(Ii2a(end)),accCI2(Ii2b(end)),accCI2(Ii2c(end));accCI3(Ii3a(end)),accCI3(Ii3b(end)),accCI3(Ii3c(end))];
TopIndicesa(:,:,expSample)=[accCI1a(Ii1aa(end)),accCI1a(Ii1ba(end)),accCI1a(Ii1ca(end));accCI2a(Ii2aa(end)),accCI2a(Ii2ba(end)),accCI2a(Ii2ca(end));accCI3a(Ii3aa(end)),accCI3a(Ii3ba(end)),accCI3a(Ii3ca(end))];
%MaxIK=[TopIK(II(1),1),TopIK(II(2),2),TopIK(II(3),3)]
[MaxIKa,IIa]=max(TopIKa);
TopIK=[TopIK1,TopIK2,TopIK3]';
[MaxIK,II]=max(TopIK);

TopSTDs=[TopSTD1,TopSTD2,TopSTD3]';
TopSTDsa=[TopSTD1a,TopSTD2a,TopSTD3a]';


MaxSTDs=[TopSTDs(II(1),1),TopSTDs(II(2),2),TopSTDs(II(3),3)];
MaxSTDsa=[TopSTDsa(II(1),1),TopSTDsa(II(2),2),TopSTDsa(II(3),3)];

AllMaxIK=[AllMaxIK;MaxIK];
AllMaxSTDs=[AllMaxSTDs;MaxSTDs];
AllMaxIKa=[AllMaxIKa;MaxIKa];
AllMaxSTDsa=[AllMaxSTDsa;MaxSTDsa];



IK2D=[IK2D,TopIK(:,1)];
IK2Da=[IK2Da,TopIKa(:,1)];
IKSTD2D=[IKSTD2D,TopSTDs(:,1)];
IKSTD2Da=[IKSTD2Da,TopSTDsa(:,1)];

end

%IK2D=[];
%IK2Da=[];
%IKSTD2D=[];
%IKSTD2Da=[];



% figure()
% bar(AllMaxIK')
% ylim([0 4.5])
% grid on
% 
% figure()
% bar(AllMaxIKa')
% ylim([0 4.5])
% grid on

figure()
%AMIK_temp=[AllMaxIK(:,1);AllMaxIK(:,2);AllMaxIK(:,3)];
%AMIK_tempa=[AllMaxIKa(:,1);AllMaxIKa(:,2);AllMaxIKa(:,3)];

IK2D=IK2D';
IK2Da=IK2Da';
IKSTD2D=IKSTD2D';
IKSTD2Da=IKSTD2Da';

AMIK_temp=[IK2D(:,1);IK2D(:,2);IK2D(:,3)];
AMIK_tempa=[IK2Da(:,1);IK2Da(:,2);IK2Da(:,3)];


IKstacked=[AMIK_temp,AMIK_tempa-AMIK_temp];
% bar(IKstacked,'stacked')
if size(AllInformation1,3)==4
IKstackedgaps=[IKstacked(1:4,:);[NaN NaN];IKstacked(5:8,:);[NaN NaN];IKstacked(9:12,:)];
else
IKstackedgaps=[IKstacked(1:3,:);[NaN NaN];IKstacked(4:6,:);[NaN NaN];IKstacked(7:9,:)];
end

IKstacked1=[IKstacked(1,:);[NaN NaN];[NaN NaN];[NaN NaN];IKstacked(4,:);[NaN NaN];[NaN NaN];[NaN NaN];IKstacked(7,:);[NaN NaN];[NaN NaN]];
IKstacked2=[[NaN NaN];IKstacked(2,:);[NaN NaN];[NaN NaN];[NaN NaN];IKstacked(5,:);[NaN NaN];[NaN NaN];[NaN NaN];IKstacked(8,:);[NaN NaN]];
IKstacked3=[[NaN NaN];[NaN NaN];IKstacked(3,:);[NaN NaN];[NaN NaN];[NaN NaN];IKstacked(6,:);[NaN NaN];[NaN NaN];[NaN NaN];IKstacked(9,:)];


figure()

hb1=bar(1:11,IKstacked1,'stacked');
hold on
hb2=bar(1:11,IKstacked2,'stacked');
hb3=bar(1:11,IKstacked3,'stacked');
colormap summer
cmap=colormap;
set(hb1,{'FaceColor'},{cmap(1,:);0.5*[1 1 1]+0.5*cmap(1,:)});
set(hb2,{'FaceColor'},{cmap(32,:);0.5*[1 1 1]+0.5*cmap(32,:)});
set(hb3,{'FaceColor'},{cmap(64,:);0.5*[1 1 1]+0.5*cmap(64,:)});

set(gcf,'Position',[1 40 560 420])
ylim([0 4.5])
%hbar=bar(IKstackedgaps,'stacked');
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.1f'))
grid on
%set(hbar,{'FaceColor'},{cmap(1,:);cmap(64,:)});
hold on
xlim([0, 12])

errorbar([-4,1:3,8],[-3;IK2D(:,1);-3],[0;IKSTD2D(:,1);0],'linestyle','none','color',.1*[1 1 1],'linewidth',2)
errorbar([0,5:7,12],[-3;IK2D(:,2);-3],[0;IKSTD2D(:,2);0],'linestyle','none','color',.1*[1 1 1],'linewidth',2)
errorbar([4,9:11,16],[-3;IK2D(:,3);-3],[0;IKSTD2D(:,3);0],'linestyle','none','color',.1*[1 1 1],'linewidth',2)

errorbar([-4,1:3,8],[-3;IK2Da(:,1);-3],[0;IKSTD2Da(:,1);0],'linestyle','none','color',.1*[1 1 1],'linewidth',2)
errorbar([0,5:7,12],[-3;IK2Da(:,2);-3],[0;IKSTD2Da(:,2);0],'linestyle','none','color',.1*[1 1 1],'linewidth',2)
errorbar([4,9:11,16],[-3;IK2Da(:,3);-3],[0;IKSTD2Da(:,3);0],'linestyle','none','color',.1*[1 1 1],'linewidth',2)

 set(gca,'xticklabel',{'24','12','6',''})
xlabel('Number of rows initialised into','interpreter','latex')
ylabel('$I_{KL}$','interpreter','latex')
text(1.5,4.7,'One SS','interpreter','latex')
text(5.5,4.7,'Two SS','interpreter','latex')
text(9.5,4.7,'Three SS','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

set(gcf,'color',[0 0 0])


% 
% if size(AllInformation1,3)==4
%     errorbar([-4,1:4,9],[-3;AllMaxIK(:,1);-3],[0;AllMaxSTDs(:,1);0],'linestyle','none','color',[0.4 0.4 0.9],'linewidth',2)
%     errorbar([1,6:9,14],[-3;AllMaxIK(:,2);-3],[0;AllMaxSTDs(:,2);0],'linestyle','none','color',[0.4 0.4 0.9],'linewidth',2)
%     errorbar([6,11:14,19],[-3;AllMaxIK(:,3);-3],[0;AllMaxSTDs(:,3);0],'linestyle','none','color',[0.4 0.4 0.9],'linewidth',2)
% else
%     errorbar([-4,1:3,8],[-3;AllMaxIK(:,1);-3],[0;AllMaxSTDs(:,1);0],'linestyle','none','color',[0.4 0.4 0.9],'linewidth',2)
%     errorbar([0,5:7,12],[-3;AllMaxIK(:,2);-3],[0;AllMaxSTDs(:,2);0],'linestyle','none','color',[0.4 0.4 0.9],'linewidth',2)
%     errorbar([4,9:11,16],[-3;AllMaxIK(:,3);-3],[0;AllMaxSTDs(:,3);0],'linestyle','none','color',[0.4 0.4 0.9],'linewidth',2)
% end
% if size(AllInformation1,3)==4
%     errorbar([-4,1:4,9],[-3;AllMaxIKa(:,1);-3],[0;AllMaxSTDsa(:,1);0],'linestyle','none','color',[0.6 0.9 0.6],'linewidth',2)
%     errorbar([1,6:9,14],[-3;AllMaxIKa(:,2);-3],[0;AllMaxSTDsa(:,2);0],'linestyle','none','color',[0.6 0.9 0.6],'linewidth',2)
%     errorbar([6,11:14,19],[-3;AllMaxIKa(:,3);-3],[0;AllMaxSTDsa(:,3);0],'linestyle','none','color',[0.6 0.9 0.6],'linewidth',2)
% else
%     errorbar([-4,1:3,8],[-3;AllMaxIKa(:,1);-3],[0;AllMaxSTDsa(:,1);0],'linestyle','none','color',[0.6 0.9 0.6],'linewidth',2)
%     errorbar([0,5:7,12],[-3;AllMaxIKa(:,2);-3],[0;AllMaxSTDsa(:,2);0],'linestyle','none','color',[0.6 0.9 0.6],'linewidth',2)
%     errorbar([4,9:11,16],[-3;AllMaxIKa(:,3);-3],[0;AllMaxSTDsa(:,3);0],'linestyle','none','color',[0.6 0.9 0.6],'linewidth',2)
% end
% ylim([0 4.5])
% xlim([0 3*(size(AllInformation1,3)+1)])
% if Experiment==1
%  %   set(gca,'xticklabel',{'24','12','6','3',''})
%     set(gca,'xticklabel',{'24','12','6',''})
%     xlabel('Number of rows initialised into','interpreter','latex')
%     ylabel('$I_{KL}$','interpreter','latex')
%     text(6.5,4.7,'$P_m$','interpreter','latex')
%     text(11.5,4.7,'$P_p$','interpreter','latex')
% elseif Experiment==2
%     set(gca,'xticklabel',{'24','48','72',''})
%     xlabel('Initial number of cells','interpreter','latex')
%     ylabel('$I_{KL}$','interpreter','latex')
%     text(5.5,4.7,'$P_m$','interpreter','latex')
%     text(9.5,4.7,'$P_p$','interpreter','latex')
% elseif Experiment==3
%     set(gca,'xticklabel',{'24','48','72',''})
%     xlabel('Initial number of cells','interpreter','latex')
%     ylabel('$I_{KL}$','interpreter','latex')
%     text(5.5,4.7,'$P_m$','interpreter','latex')
%     text(9.5,4.7,'$P_p$','interpreter','latex')
% end
% 
% 
% % 
% % figure()
% % IKstacked2(:,:,1)=reshape(IKstacked(:,1),[4,3])';
% % IKstacked2(:,:,2)=reshape(IKstacked(:,2),[4,3])';
% % plotBarStackGroups(IKstacked2,{'1','2','3'})
% end
% for fignum=1:3
%     figure(fignum)
%     pause(1)
% %    export_fig( gcf, ['Experiment',num2str(fignum)], '-painters','-nocrop','-transparent','-pdf' );
% end
    
