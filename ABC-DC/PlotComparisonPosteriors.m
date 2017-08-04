%% Take the ABC DC chains and the ABC rejection theta to obtain posteriors. 
%% Saving out the plots using print function is at bottom of this file

%% Compute Densities:

load('ABCsimulationsCXYSD6DF1.mat');
load('ABCexperimentalCXYSD6DF1.mat');

%% ABC-DC Density:
thetaRecC=thetaRecord;

Pm_xi=linspace(0,1,1001);
Pp_xi=linspace(0,0.01,1001);

Pm_Fd(:,1)=ksdensity(thetaRecC(1,4e3+1:5e3),Pm_xi);
Pp_Fd(:,1)=ksdensity(thetaRecC(2,4e3+1:5e3),Pp_xi);

Pm_Fd(:,2)=ksdensity(thetaRecC(1,5e3+1:6e3),Pm_xi);
Pp_Fd(:,2)=ksdensity(thetaRecC(2,5e3+1:6e3),Pp_xi);

Pm_Fd(:,3)=ksdensity(thetaRecC(1,6e3+1:7e3),Pm_xi);
Pp_Fd(:,3)=ksdensity(thetaRecC(2,6e3+1:7e3),Pp_xi);

Pm_Fd(:,4)=ksdensity(thetaRecC(1,7e3+1:8e3),Pm_xi);
Pp_Fd(:,4)=ksdensity(thetaRecC(2,7e3+1:8e3),Pp_xi);

Pm_Fd(:,5)=ksdensity(thetaRecC(1,8e3+1:10e3),Pm_xi);
Pp_Fd(:,5)=ksdensity(thetaRecC(2,8e3+1:10e3),Pp_xi);

%% ABC-rej Density:

k=200;
thetanew1=theta(:,I2c(1:k));
thetanew2=RegressionAdjustmentFunction(theta(:,I2c(1:k)),Ssim(1:24,I2c(1:k)),Sexp(1:24),Ssigma(1:24),0);
thetanew3=RegressionAdjustmentFunction(theta(:,I2c(1:k)),Ssim(1:24,I2c(1:k)),Sexp(1:24),Ssigma(1:24),0.25);
thetanew4=RegressionAdjustmentFunction(theta(:,I2c(1:k)),Ssim(1:24,I2c(1:k)),Sexp(1:24),Ssigma(1:24),0.2);

for i=1:4
    switch i
        case 1
            thetan=thetanew1;
        case 2
            thetan=thetanew2;
        case 3
            thetan=thetanew3;
        case 4
            thetan=thetanew4;
    end
%figure()
%subplot(2,1,1)
%plot(theta(1,I2c(1:k)))
%hold on, plot(thetan(1,1:k))
%figure()
%subplot(2,1,2)
%plot(theta(2,I2c(1:k)))
%hold on, plot(thetan(2,1:k))

% ks density
Pm_xi=linspace(0,1,1001);
Pp_xi=linspace(0,0.01,1001);

Pm_Fr(:,i)=ksdensity(thetan(1,:),Pm_xi);
Pp_Fr(:,i)=ksdensity(thetan(2,:),Pp_xi);
end

%% Now perform plotting of selected densities:

cmapshort=get(gca,'colororder');

% 1st is the regression adjusted ABC rejection:
Pm_F(:,1)=Pm_Fr(:,3); 
Pp_F(:,1)=Pp_Fr(:,3); 
% 2nd is the regression adjusted ABC MCMC:
Pm_F(:,2)=Pm_Fd(:,2); 
Pp_F(:,2)=Pp_Fd(:,2); 
% 3rd is the ABC-DC K=4:
Pm_F(:,3)=Pm_Fd(:,5); 
Pp_F(:,3)=Pp_Fd(:,5); 


figure()
%subplot(2,1,1)
linecol=[1 1 1]*0;
shadedplot(Pm_xi,Pm_F(:,3)',zeros(1001,1)',cmapshort(3,:),linecol)
hold on
shadedplot(Pm_xi,Pm_F(:,2)',zeros(1001,1)',cmapshort(2,:),linecol)
hold on
shadedplot(Pm_xi,Pm_F(:,1)',zeros(1001,1)',cmapshort(1,:),linecol)
alpha(.7)
xlim([0 0.5])
xlabel('$P_m$','interpreter','latex')
ylabel('Posterior','interpreter','latex')
prettymyplots2(6,4,'%0.1f','%2.0f')

%subplot(2,1,2)
figure()
linecol=[1 1 1]*0;
shadedplot(Pp_xi,Pp_F(:,3)',zeros(1001,1)',cmapshort(3,:),linecol)
hold on
shadedplot(Pp_xi,Pp_F(:,2)',zeros(1001,1)',cmapshort(2,:),linecol)
hold on
shadedplot(Pp_xi,Pp_F(:,1)',zeros(1001,1)',cmapshort(1,:),linecol)
alpha(.7)
xlim([0.15 0.35]/100)
ylim([0 8000])
xlabel('$P_p$','interpreter','latex')
ylabel('Posterior','interpreter','latex')
prettymyplots2(5,5,'%1.1f','%1.1f')

%% Save plots using the below

% h=figure(4)
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% print(h,'filename','-dpdf','-r0')
