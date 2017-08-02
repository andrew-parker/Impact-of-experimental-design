% %% ABC rejection quantiles:
% % ABC rejection spits out a theta (2x20000) and an I (1x20000)
% % This means that theta(1,I) is the sorted distances
% % So as we only accept the smallest 1%, we can find the 90% CI using:
% load('matsABCrejSD24DF1.mat')
% % qrN24 is quantiles for ABCrej with 24 rows for ss N:
% qrN24(1,:)=quantile(thN(1,IN(1:200)),[0.1 0.5 0.9]);
% qrN24(2,:)=quantile(thN(2,IN(1:200)),[0.1 0.5 0.9]);
% qrX24(1,:)=quantile(thX(1,IX(1:200)),[0.1 0.5 0.9]);
% qrX24(2,:)=quantile(thX(2,IX(1:200)),[0.1 0.5 0.9]);
% qrC24(1,:)=quantile(thC(1,IC(1:200)),[0.1 0.5 0.9]);
% qrC24(2,:)=quantile(thC(2,IC(1:200)),[0.1 0.5 0.9]);
% % now we do the same again for 6 rows:
% load('matsABCrejSD6DF1.mat')
% qrN6(1,:)=quantile(thN(1,IN(1:200)),[0.1 0.5 0.9]);
% qrN6(2,:)=quantile(thN(2,IN(1:200)),[0.1 0.5 0.9]);
% qrX6(1,:)=quantile(thX(1,IX(1:200)),[0.1 0.5 0.9]);
% qrX6(2,:)=quantile(thX(2,IX(1:200)),[0.1 0.5 0.9]);
% qrC6(1,:)=quantile(thC(1,IC(1:200)),[0.1 0.5 0.9]);
% qrC6(2,:)=quantile(thC(2,IC(1:200)),[0.1 0.5 0.9]);
% qrN6(2,:)=qrN6(2,:)*100;
% qrN24(2,:)=qrN24(2,:)*100;
% qrC6(2,:)=qrC6(2,:)*100;
% qrC24(2,:)=qrC24(2,:)*100;
% qrX6(2,:)=qrX6(2,:)*100;
% qrX24(2,:)=qrX24(2,:)*100;
% 
% %% Print resulting quantiles:
% str{1,1}=([num2str(qrN6(1,2),3),' (',num2str(qrN6(1,1),3),', ',num2str(qrN6(1,3),3),')']);
% str{1,2}=([num2str(qrX6(1,2),3),' (',num2str(qrX6(1,1),3),', ',num2str(qrX6(1,3),3),')']);
% str{1,3}=([num2str(qrC6(1,2),3),' (',num2str(qrC6(1,1),3),', ',num2str(qrC6(1,3),3),')']);
% 
% str{1,4}=([num2str(qrN24(1,2),3),' (',num2str(qrN24(1,1),3),', ',num2str(qrN24(1,3),3),')']);
% str{1,5}=([num2str(qrX24(1,2),3),' (',num2str(qrX24(1,1),3),', ',num2str(qrX24(1,3),3),')']);
% str{1,6}=([num2str(qrC24(1,2),3),' (',num2str(qrC24(1,1),3),', ',num2str(qrC24(1,3),3),')']);
% 
% str{3,1}=([num2str(qrN6(2,2),3),' (',num2str(qrN6(2,1),3),', ',num2str(qrN6(2,3),3),')']);
% str{3,2}=([num2str(qrX6(2,2),3),' (',num2str(qrX6(2,1),3),', ',num2str(qrX6(2,3),3),')']);
% str{3,3}=([num2str(qrC6(2,2),3),' (',num2str(qrC6(2,1),3),', ',num2str(qrC6(2,3),3),')']);
% 
% str{3,4}=([num2str(qrN24(2,2),3),' (',num2str(qrN24(2,1),3),', ',num2str(qrN24(2,3),3),')']);
% str{3,5}=([num2str(qrX24(2,2),3),' (',num2str(qrX24(2,1),3),', ',num2str(qrX24(2,3),3),')']);
% str{3,6}=([num2str(qrC24(2,2),3),' (',num2str(qrC24(2,1),3),', ',num2str(qrC24(2,3),3),')']);
% %str=str'
% %inp.tableColLabels = {'$P_m$','$P_p$'};
% %inp.tableRowLabels = {'N 6','X 6','C_{XY} 6','N 24','X 24','C_{XY} 24'};
% %inp.data=str;
% %latexTable(inp)
% %% ABC DC quantiles:
% % This depends on for which population we want to compute the quantile
% % The scheme is decreasing delta for 6 populations, then 1 population with regression adjustment applied, then increase K 3 times
% % So the final population with 4 clones is thetaRecord(1,8001:10000); But
% % we can also compare with 1 clone and the smallest delta: (ABC MCMC) thetaRecord(1,6001:7000)
% load('matsABCDCSD24DF1.mat')
% qdN24(1,:)=quantile(thetaRecN(1,8001:10000),[0.1 0.5 0.9]);
% qdN24(2,:)=quantile(thetaRecN(2,8001:10000),[0.1 0.5 0.9]);
% qdX24(1,:)=quantile(thetaRecX(1,8001:10000),[0.1 0.5 0.9]);
% qdX24(2,:)=quantile(thetaRecX(2,8001:10000),[0.1 0.5 0.9]);
% qdC24(1,:)=quantile(thetaRecC(1,8001:10000),[0.1 0.5 0.9]);
% qdC24(2,:)=quantile(thetaRecC(2,8001:10000),[0.1 0.5 0.9]);
% load('matsABCDCSD6DF1.mat')
% qdN6(1,:)=quantile(thetaRecN(1,8001:10000),[0.1 0.5 0.9]);
% qdN6(2,:)=quantile(thetaRecN(2,8001:10000),[0.1 0.5 0.9]);
% qdX6(1,:)=quantile(thetaRecX(1,8001:10000),[0.1 0.5 0.9]);
% qdX6(2,:)=quantile(thetaRecX(2,8001:10000),[0.1 0.5 0.9]);
% qdC6(1,:)=quantile(thetaRecC(1,8001:10000),[0.1 0.5 0.9]);
% qdC6(2,:)=quantile(thetaRecC(2,8001:10000),[0.1 0.5 0.9]);
% qdN6(2,:)=qdN6(2,:)*100;
% qdN24(2,:)=qdN24(2,:)*100;
% qdC6(2,:)=qdC6(2,:)*100;
% qdC24(2,:)=qdC24(2,:)*100;
% qdX6(2,:)=qdX6(2,:)*100;
% qdX24(2,:)=qdX24(2,:)*100;
% 
% %% Print resulting quantiles:
% %clear str
% 
% str{2,1}=([num2str(qdN6(1,2),3),' (',num2str(qdN6(1,1),3),', ',num2str(qdN6(1,3),3),')']);
% str{2,2}=([num2str(qdX6(1,2),3),' (',num2str(qdX6(1,1),3),', ',num2str(qdX6(1,3),3),')']);
% str{2,3}=([num2str(qdC6(1,2),3),' (',num2str(qdC6(1,1),3),', ',num2str(qdC6(1,3),3),')']);
% 
% str{2,4}=([num2str(qdN24(1,2),3),' (',num2str(qdN24(1,1),3),', ',num2str(qdN24(1,3),3),')']);
% str{2,5}=([num2str(qdX24(1,2),3),' (',num2str(qdX24(1,1),3),', ',num2str(qdX24(1,3),3),')']);
% str{2,6}=([num2str(qdC24(1,2),3),' (',num2str(qdC24(1,1),3),', ',num2str(qdC24(1,3),3),')']);
% 
% str{4,1}=([num2str(qdN6(2,2),3),' (',num2str(qdN6(2,1),3),', ',num2str(qdN6(2,3),3),')']);
% str{4,2}=([num2str(qdX6(2,2),3),' (',num2str(qdX6(2,1),3),', ',num2str(qdX6(2,3),3),')']);
% str{4,3}=([num2str(qdC6(2,2),3),' (',num2str(qdC6(2,1),3),', ',num2str(qdC6(2,3),3),')']);
% 
% str{4,4}=([num2str(qdN24(2,2),3),' (',num2str(qdN24(2,1),3),', ',num2str(qdN24(2,3),3),')']);
% str{4,5}=([num2str(qdX24(2,2),3),' (',num2str(qdX24(2,1),3),', ',num2str(qdX24(2,3),3),')']);
% str{4,6}=([num2str(qdC24(2,2),3),' (',num2str(qdC24(2,1),3),', ',num2str(qdC24(2,3),3),')']);
% str=str'
% % 6 rows, 2 columns
% inp.tableColLabels = {'$P_m$ (ABC-rej)','$P_m$ (ABC-DC)','$P_p$ (ABC-DC)','$P_p$ (ABC-DC)'};
% inp.tableRowLabels = {'$N$ 6','$X$ 6','$C_{XY}$ 6','$N$ 24','$X$ 24','$C_{XY}$ 24'};
% inp.data=str;
% latexTable(inp)



%%
%%
%% Print maximum likelihood estimates of posteriors (using ksdensity): 
Pm_xi=linspace(0,0.99,10001);
Pp_xi=linspace(0,0.01,10001);
k=200;

%% N Pp
load('matsABCrejNSD6DF1.mat')
load('matsABCDCNSD6DF1.mat')
L=1;
I=I1c;
theta=RegressionAdjustmentFunction(theta(:,I(1:k)),Ssim(1:L,I(1:k)),Sexp(1:L),Ssigma(1:L),0);
I=1:k;
Pp_KDE=ksdensity(theta(2,I(1:k)),Pp_xi,'support',[0 0.01]);
Pp_MLE=100*Pp_xi(Pp_KDE==max(Pp_KDE));
Pp_QTL=100*quantile(theta(2,:),[0.05,0.95]);
Pp_DCMLE=100*mean(thetaRecord(2,8001:10000));
Pp_DCQTL=100*quantile(thetaRecord(2,8001:10000),[0.05,0.95]);

str{4,1}=([num2str(Pp_MLE,3),' (',num2str(Pp_QTL(1),3),', ',num2str(Pp_QTL(2),3),')']);
str{4,2}=([num2str(Pp_DCMLE,3),' (',num2str(Pp_DCQTL(1),3),', ',num2str(Pp_DCQTL(2),3),')']);
str{4,3}=([num2str((Tsim+Tcal)/(TK+Ttol),3)]);

%% CXY Pm
load('matsABCrejCXYSD6DF1.mat')
load('matsABCDCCXYSD6DF1.mat')
L=24;
I=I2c;
theta=RegressionAdjustmentFunction(theta(:,I(1:k)),Ssim(1:L,I(1:k)),Sexp(1:L),Ssigma(1:L),0);
I=1:k;
%densityplot(theta2,1:k,'RA')
Pm_KDE=ksdensity(theta(1,I(1:k)),Pm_xi,'support',[0 0.99]);
Pm_MLE=Pm_xi(Pm_KDE==max(Pm_KDE));
Pm_QTL=quantile(theta(1,:),[0.05,0.95]);
Pm_DCMLE=mean(thetaRecord(1,8001:10000));
Pm_DCQTL=quantile(thetaRecord(1,8001:10000),[0.05,0.95]);

Pp_KDE=ksdensity(theta(2,I(1:k)),Pp_xi,'support',[0 0.01]);
Pp_MLE=100*Pp_xi(Pp_KDE==max(Pp_KDE));
Pp_QTL=100*quantile(theta(2,:),[0.05,0.95]);
Pp_DCMLE=100*mean(thetaRecord(2,8001:10000));
Pp_DCQTL=100*quantile(thetaRecord(2,8001:10000),[0.05,0.95]);

str{3,1}=([num2str(Pp_MLE,3),' (',num2str(Pp_QTL(1),3),', ',num2str(Pp_QTL(2),3),')']);
str{3,2}=([num2str(Pp_DCMLE,3),' (',num2str(Pp_DCQTL(1),3),', ',num2str(Pp_DCQTL(2),3),')']);
str{3,3}=([num2str((Tsim+Tcal)/(TK+Ttol),3)]);

str{2,1}=([num2str(Pm_MLE,3),' (',num2str(Pm_QTL(1),3),', ',num2str(Pm_QTL(2),3),')']);
str{2,2}=([num2str(Pm_DCMLE,3),' (',num2str(Pm_DCQTL(1),3),', ',num2str(Pm_DCQTL(2),3),')']);
str{2,3}=([num2str((Tsim+Tcal)/(TK+Ttol),3)]);


load('matsABCrejXSD6DF1.mat')
load('matsABCDCXSD6DF1.mat')
L=1;
I=I6A;

theta=RegressionAdjustmentFunction(theta(:,I(1:k)),Ssim(1:L,I(1:k)),Sexp(1:L),Ssigma(1:L),0);
I=1:k;
%densityplot(theta2,1:k,'RA')
Pm_KDE=ksdensity(theta(1,I(1:k)),Pm_xi,'support',[0 0.99]);
Pm_MLE=Pm_xi(Pm_KDE==max(Pm_KDE));
Pm_QTL=quantile(theta(1,:),[0.05,0.95]);
Pm_DCMLE=mean(thetaRecord(1,8001:10000));
Pm_DCQTL=quantile(thetaRecord(1,8001:10000),[0.05,0.95]);

str{1,1}=([num2str(Pm_MLE,3),' (',num2str(Pm_QTL(1),3),', ',num2str(Pm_QTL(2),3),')']);
str{1,2}=([num2str(Pm_DCMLE,3),' (',num2str(Pm_DCQTL(1),3),', ',num2str(Pm_DCQTL(2),3),')']);
str{1,3}=([num2str((Tsim+Tcal)/(TK+Ttol),3)]);

inp.tableColLabels = {'ABC-rej MLE (90\% CI)','ABC-DC MLE (90\% CI)','Speedup'};
inp.tableRowLabels = {'$\|x\|$: $P_m$','$C_{XY}$: $P_m$','$C_{XY}$: $P_p(\times 10^{-2})$','$N$: $P_p(\times 10^{-2})$'};
inp.data=str;
latexTable(inp)


%%
% %%
% %% Print Means and Covariances instead
% clear all
% 
% load('matsABCrejSD24DF1.mat')
% % qrN24 is quantiles for ABCrej with 24 rows for ss N:
% qrN24(1,:)=[mean(thN(1,IN(1:200))),std(thN(1,IN(1:200)))]';
% qrN24(2,:)=[mean(thN(2,IN(1:200))),std(thN(2,IN(1:200)))];
% qrX24(1,:)=[mean(thX(1,IX(1:200))),std(thX(1,IX(1:200)))];
% qrX24(2,:)=[mean(thX(2,IX(1:200))),std(thX(2,IX(1:200)))];
% qrC24(1,:)=[mean(thC(1,IC(1:200))),std(thC(1,IC(1:200)))];
% qrC24(2,:)=[mean(thC(2,IC(1:200))),std(thC(2,IC(1:200)))];
% % now we do the same again for 6 rows:
% load('matsABCrejSD6DF1.mat')
% qrN6(1,:)=[mean(thN(1,IN(1:200))),std(thN(1,IN(1:200)))];
% qrN6(2,:)=[mean(thN(2,IN(1:200))),std(thN(2,IN(1:200)))];
% qrX6(1,:)=[mean(thX(1,IX(1:200))),std(thX(1,IX(1:200)))];
% qrX6(2,:)=[mean(thX(2,IX(1:200))),std(thX(2,IX(1:200)))];
% qrC6(1,:)=[mean(thC(1,IC(1:200))),std(thC(1,IC(1:200)))];
% qrC6(2,:)=[mean(thC(2,IC(1:200))),std(thC(2,IC(1:200)))];
% qrN6(2,:)=qrN6(2,:)*100;
% qrN24(2,:)=qrN24(2,:)*100;
% qrC6(2,:)=qrC6(2,:)*100;
% qrC24(2,:)=qrC24(2,:)*100;
% qrX6(2,:)=qrX6(2,:)*100;
% qrX24(2,:)=qrX24(2,:)*100;
% 
% str{1,1}=([num2str(qrN6(1,1),3),' (',num2str(qrN6(1,2),3),')']);
% str{1,2}=([num2str(qrX6(1,1),3),' (',num2str(qrX6(1,2),3),')']);
% str{1,3}=([num2str(qrC6(1,1),3),' (',num2str(qrC6(1,2),3),')']);
% str{1,4}=([num2str(qrN24(1,1),3),' (',num2str(qrN24(1,2),3),')']);
% str{1,5}=([num2str(qrX24(1,1),3),' (',num2str(qrX24(1,2),3),')']);
% str{1,6}=([num2str(qrC24(1,1),3),' (',num2str(qrC24(1,2),3),')']);
% str{3,1}=([num2str(qrN6(2,1),3),' (',num2str(qrN6(2,2),3),')']);
% str{3,2}=([num2str(qrX6(2,1),3),' (',num2str(qrX6(2,2),3),')']);
% str{3,3}=([num2str(qrC6(2,1),3),' (',num2str(qrC6(2,2),3),')']);
% str{3,4}=([num2str(qrN24(2,1),3),' (',num2str(qrN24(2,2),3),')']);
% str{3,5}=([num2str(qrX24(2,1),3),' (',num2str(qrX24(2,2),3),')']);
% str{3,6}=([num2str(qrC24(2,1),3),' (',num2str(qrC24(2,2),3),')']);
% 
% load('matsABCDCSD24DF1.mat')
% qdN24(1,:)=[mean(thetaRecN(1,8001:10000)),2*std(thetaRecN(1,8001:10000))];
% qdN24(2,:)=[mean(thetaRecN(2,8001:10000)),2*std(thetaRecN(2,8001:10000))];
% qdX24(1,:)=[mean(thetaRecX(1,8001:10000)),2*std(thetaRecX(1,8001:10000))];
% qdX24(2,:)=[mean(thetaRecX(2,8001:10000)),2*std(thetaRecX(2,8001:10000))];
% qdC24(1,:)=[mean(thetaRecC(1,8001:10000)),2*std(thetaRecC(1,8001:10000))];
% qdC24(2,:)=[mean(thetaRecC(2,8001:10000)),2*std(thetaRecC(2,8001:10000))];
% load('matsABCDCSD6DF1.mat')
% qdN6(1,:)=[mean(thetaRecN(1,8001:10000)),2*std(thetaRecN(1,8001:10000))];
% qdN6(2,:)=[mean(thetaRecN(2,8001:10000)),2*std(thetaRecN(2,8001:10000))];
% qdX6(1,:)=[mean(thetaRecX(1,8001:10000)),2*std(thetaRecX(1,8001:10000))];
% qdX6(2,:)=[mean(thetaRecX(2,8001:10000)),2*std(thetaRecX(2,8001:10000))];
% qdC6(1,:)=[mean(thetaRecC(1,8001:10000)),2*std(thetaRecC(1,8001:10000))];
% qdC6(2,:)=[mean(thetaRecC(2,8001:10000)),2*std(thetaRecC(2,8001:10000))];
% qdN6(2,:)=qdN6(2,:)*100;
% qdN24(2,:)=qdN24(2,:)*100;
% qdC6(2,:)=qdC6(2,:)*100;
% qdC24(2,:)=qdC24(2,:)*100;
% qdX6(2,:)=qdX6(2,:)*100;
% qdX24(2,:)=qdX24(2,:)*100;
% 
% str{2,1}=([num2str(qdN6(1,1),3),' (',num2str(qdN6(1,2),3),')']);
% str{2,2}=([num2str(qdX6(1,1),3),' (',num2str(qdX6(1,2),3),')']);
% str{2,3}=([num2str(qdC6(1,1),3),' (',num2str(qdC6(1,2),3),')']);
% str{2,4}=([num2str(qdN24(1,1),3),' (',num2str(qdN24(1,2),3),')']);
% str{2,5}=([num2str(qdX24(1,1),3),' (',num2str(qdX24(1,2),3),')']);
% str{2,6}=([num2str(qdC24(1,1),3),' (',num2str(qdC24(1,2),3),')']);
% str{4,1}=([num2str(qdN6(2,1),3),' (',num2str(qdN6(2,2),3),')']);
% str{4,2}=([num2str(qdX6(2,1),3),' (',num2str(qdX6(2,2),3),')']);
% str{4,3}=([num2str(qdC6(2,1),3),' (',num2str(qdC6(2,2),3),')']);
% str{4,4}=([num2str(qdN24(2,1),3),' (',num2str(qdN24(2,2),3),')']);
% str{4,5}=([num2str(qdX24(2,1),3),' (',num2str(qdX24(2,2),3),')']);
% str{4,6}=([num2str(qdC24(2,1),3),' (',num2str(qdC24(2,2),3),')']);
% str=str'
% 
% %inp.tableColLabels = {'$P_m$ (ABC-rej)','$P_m$ (ABC-DC)','$P_p$ (ABC-rej)','$P_p$ (ABC-DC)'};
% inp.tableColLabels = {'ABC-rej','ABC-DC','ABC-rej','ABC-DC'};
% inp.tableRowLabels = {'$N$ 6','$X$ 6','$C_{XY}$ 6','$N$ 24','$X$ 24','$C_{XY}$ 24'};
% inp.data=str;
% latexTable(inp)
% 
% 
% %%
% %% Bootstrap
% 
% 
% Pmmin=0;
% Pmmax=0.99;
% datapoints=1000;
% xi=linspace(Pmmin,Pmmax,datapoints);
% for j=1:8
%     dataset=thetaRecN(1,(j-1)*1000+1:j*1000);
%     pdfdata(j,:)=ksdensity(dataset,xi,'support',[Pmmin Pmmax]);
% end
% dataset=thetaRecN(1,8001:10000);
% pdfdata(9,:)=ksdensity(dataset,xi,'support',[Pmmin Pmmax]);
% 
% 
% 
% m=bootstrp(100,@mean,dataset);
% 
% 
% % %% ABC MCMC quantiles
% % load('matsABCDCSD24DF1.mat')
% % qdN24(1,:)=quantile(thetaRecN(1,5001:6000),[0.1 0.5 0.9]);
% % qdN24(2,:)=quantile(thetaRecN(2,5001:6000),[0.1 0.5 0.9]);
% % qdX24(1,:)=quantile(thetaRecX(1,5001:6000),[0.1 0.5 0.9]);
% % qdX24(2,:)=quantile(thetaRecX(2,5001:6000),[0.1 0.5 0.9]);
% % qdC24(1,:)=quantile(thetaRecC(1,5001:6000),[0.1 0.5 0.9]);
% % qdC24(2,:)=quantile(thetaRecC(2,5001:6000),[0.1 0.5 0.9]);
% % load('matsABCDCSD6DF1.mat')
% % qdN6(1,:)=quantile(thetaRecN(1,5001:6000),[0.1 0.5 0.9]);
% % qdN6(2,:)=quantile(thetaRecN(2,5001:6000),[0.1 0.5 0.9]);
% % qdX6(1,:)=quantile(thetaRecX(1,5001:6000),[0.1 0.5 0.9]);
% % qdX6(2,:)=quantile(thetaRecX(2,5001:6000),[0.1 0.5 0.9]);
% % qdC6(1,:)=quantile(thetaRecC(1,5001:6000),[0.1 0.5 0.9]);
% % qdC6(2,:)=quantile(thetaRecC(2,5001:6000),[0.1 0.5 0.9]);
% 
% %% Plot:
% figure()
% errorbar(qdN6(1,2),qdN6(2,2),qdN6(2,2)-qdN6(2,1),qdN6(2,3)-qdN6(2,2),qdN6(1,2)-qdN6(1,1),qdN6(1,3)-qdN6(1,2))
% hold on
% errorbar(qdX6(1,2),qdX6(2,2),qdX6(2,2)-qdX6(2,1),qdX6(2,3)-qdX6(2,2),qdX6(1,2)-qdX6(1,1),qdX6(1,3)-qdX6(1,2))
% errorbar(qdC6(1,2),qdC6(2,2),qdC6(2,2)-qdC6(2,1),qdC6(2,3)-qdC6(2,2),qdC6(1,2)-qdC6(1,1),qdC6(1,3)-qdC6(1,2))
%     xlim([0 1])
%     ylim([0 0.01])
% errorbar(qdN24(1,2),qdN24(2,2),qdN24(2,2)-qdN24(2,1),qdN24(2,3)-qdN24(2,2),qdN24(1,2)-qdN24(1,1),qdN24(1,3)-qdN24(1,2))
% %hold on
% errorbar(qdX24(1,2),qdX24(2,2),qdX24(2,2)-qdX24(2,1),qdX24(2,3)-qdX24(2,2),qdX24(1,2)-qdX24(1,1),qdX24(1,3)-qdX24(1,2))
% errorbar(qdC24(1,2),qdC24(2,2),qdC24(2,2)-qdC24(2,1),qdC24(2,3)-qdC24(2,2),qdC24(1,2)-qdC24(1,1),qdC24(1,3)-qdC24(1,2))
% 
% %% Print resulting quantiles:
% disp([num2str(qdN6(1,2),3),' (',num2str(qdN6(1,1),3),', ',num2str(qdN6(1,3),3),')'])
% disp([num2str(qdX6(1,2),3),' (',num2str(qdX6(1,1),3),', ',num2str(qdX6(1,3),3),')'])
% disp([num2str(qdC6(1,2),3),' (',num2str(qdC6(1,1),3),', ',num2str(qdC6(1,3),3),')'])
% 
% disp([num2str(qdN24(1,2),3),' (',num2str(qdN24(1,1),3),', ',num2str(qdN24(1,3),3),')'])
% disp([num2str(qdX24(1,2),3),' (',num2str(qdX24(1,1),3),', ',num2str(qdX24(1,3),3),')'])
% disp([num2str(qdC24(1,2),3),' (',num2str(qdC24(1,1),3),', ',num2str(qdC24(1,3),3),')'])
% 
% disp([num2str(qdN6(2,2),3),' (',num2str(qdN6(2,1),3),', ',num2str(qdN6(2,3),3),')'])
% disp([num2str(qdX6(2,2),3),' (',num2str(qdX6(2,1),3),', ',num2str(qdX6(2,3),3),')'])
% disp([num2str(qdC6(2,2),3),' (',num2str(qdC6(2,1),3),', ',num2str(qdC6(2,3),3),')'])
% 
% disp([num2str(qdN24(2,2),3),' (',num2str(qdN24(2,1),3),', ',num2str(qdN24(2,3),3),')'])
% disp([num2str(qdX24(2,2),3),' (',num2str(qdX24(2,1),3),', ',num2str(qdX24(2,3),3),')'])
% disp([num2str(qdC24(2,2),3),' (',num2str(qdC24(2,1),3),', ',num2str(qdC24(2,3),3),')'])
% 








