function [thetaRecord,acceptedAll]=ABC_DC_Cxy(SD,DF)
% This code was originally for doing ABC with varying scratch amounts (as
% in my poster)
% I've got rid of all mention of "replicates" too
% ABC data cloning works in two steps, 1st it reduces tolerance 2nd it increases #clones
%% Setup ABC-DC: We define a tolerance scheme and a clones scheme
%trimmedThetaExperiments=zeros(2,100,12);
replicates=10;
ROWS=24;
%deltaScheme=[10,5,2,1];
deltaScheme=[5,2,1,0.5,0.25];
P=length(deltaScheme); % Stage 1 is run P times,
clonesScheme=[1,2,3,4];
Q=length(clonesScheme); % Stage 2 is run Q times,
deltaSamplesScheme=repmat(1000,[1,P]);
%deltaSamplesScheme(end)=5000;
clonesSamplesScheme=repmat(1000,[1,Q]);
clonesSamplesScheme(Q)=2000;
%if P+Q~=length(samplesScheme)
%error('incorrect samplesScheme length')
%end
tic

load(['ABCexperimentalCXYSD',int2str(SD),'DF',int2str(DF),'.mat'])
load(['sigma_SD',num2str(SD),'DF',num2str(DF),'.mat'])
Ssigma=sigmalarge(2:25);

%% Setup Model variables and parameters:

%% 1: Compute synthetic experiments: data and statistics
% Pm=0.25;
% Pp=0.0025;
% 
% %replicates=10;
% 
% ROWS=24;
% 
% Mexp=zeros(24,32,replicates);
% for reps=1:replicates
% M=zeros(32,24);
% M(randsample(32*ScratchDensity,DensityFraction*24))=1;
% %Mexp(:,:)=M';
% Mexp(:,:,reps)=M';
% end

% Ssigma=[1 1]; 
%load(['sigma_SD',num2str(SD),'DF',num2str(DF),'.mat'])
%Ssigma=sigmalarge(2:25);

% Aexp=false(ROWS,ROWS*4/3,replicates); % Final matrices cropped
% 
% 
% %% Simulate to find experimental SS
% 
% %S1exp=0;
% %S2exp=0;
% Sexp=zeros(24,1);
% Sexpall=zeros(24,replicates);
% for reps=1:replicates
%     [Aexp(:,:,reps),~,~]=ABCmex(Pm,Pp,Mexp(:,:,reps));
%     [I,J]=find(Aexp(:,:,reps));
%     D=mandist([I,J]'); %manhattan distance, though could use boxdist.
%     for rr=1:24
%         Sexpall(rr,reps)=length(find(triu(D)==rr));
%     end
% end
% Sexp=mean(Sexpall,2);

%% Stage 1: ABC MCMC
% Initialise:
acceptedP=zeros(1,P);
acceptedQ=zeros(1,Q);
K=clonesScheme(1); % K=1;
Delta=deltaScheme(1); % Delta is large initially

% Initialise A and Aproposed with replicates #entries (as opposed to
% clones, because K=1 initially!)
A=false(ROWS,ROWS*4/3,replicates); % Final matrices
Aproposed=false(ROWS,ROWS*4/3,replicates); % Final matrices
% 
% N0=5;
% FR=zeros(N0,36,K); % One trajectory with 36 points (the cell with index 0)
% FC=zeros(N0,36,K); % One trajectory with 36 points (the cell with index 0)
% FRproposed=zeros(N0,36,K); % One trajectory with 36 points (the cell with index 0)
% FCproposed=zeros(N0,36,K); % One trajectory with 36 points (the cell with index 0)

thetaOld(1,1)=0.5; % either propose or generate from prior
thetaOld(2,1)=0.003;
thetaRecord=[];
Sregression=zeros(24,deltaSamplesScheme(P));

tic;
for p=1:P % perform ABC MCMC for every delta
    
    
    if p==1
        theta=zeros(2,deltaSamplesScheme(p));    
        theta(:,1)=thetaOld;
    else
        thetaTemp=theta(:,p-1);
        theta=zeros(2,deltaSamplesScheme(p));    
        theta(:,1)=thetaTemp;
        Delta=deltaScheme(p);
    end
    
    % Initialise SS
    % K=1, so replace K with reps here:
    Sall=zeros(24,replicates);
    Sproposedall=zeros(24,replicates);
    
    S=zeros(24,1);
    Sproposed=zeros(24,1);
    
    
    %F1=zeros(N0,1);
    %S2=zeros(1,K);
    %S2proposed=zeros(1,K);
    
    
    % Generate Data from thetaStar and compute qstar
    for reps=1:replicates
       
        % In this part of the code, K=1, so no clones
        
           % [A(:,:,i),FR(:,:,reps),FC(:,:,reps)]=ABCmex(thetaOld(1,1),thetaOld(2,1),Mexp(:,:,reps));
            [A(:,:,reps),~,~]=ABCmex(thetaOld(1,1),thetaOld(2,1),Mexp(:,:,reps));
            % Compute SS: 1 #cells
    %        S1(i)=(sum(sum(A(:,:,i))));
            % Compute SS: 2 Trajectory distance
    %        for n=1:N0
    %            Trajectory=[FR(n,:,i);FC(n,:,i)];
    %            F1(n,1)=sum(sum(abs(diff(Trajectory')')));
    %        end
    %        meanF1=mean(F1(1:N0,:),1); % average over the tracks, N0
    %        SS6A=mean(meanF1); % average over the replicates
    %        S2(i)=SS6A;
            % % Compute SS: 3 Pairwise correlation:

             [I,J]=find(A(:,:,reps));
             D=mandist([I,J]'); 
             for rr=1:24
                 Sall(rr,reps)=length(find(triu(D)==rr));
             end
    end
    S=mean(Sall,2);
    % 
    %S=[S1;S2];
    %S=S1';
    %Sexp=[S1exp;S2exp];
    %Sexp=S1exp;
    %Weights=diag([Ssigma(1),sigma(2)]);
    Weights=diag(Ssigma)^(-1);
    qstar=prod(exp(-diag((S-Sexp)'*Weights*(S-Sexp))./(24*2*Delta^2)));
    %qstar=prod(exp(-(S1-S1exp).^2/(2*Delta^2)));
    %qstar=prod(exp(-(B-Bexp).^2/(2*Delta^2)));
    
    Sregression(:,1)=S;
    for samples=1:(deltaSamplesScheme(p)-1)
        
        % Generate thetaProposed from thetaStar+Noise and compute qNew

        % thetaproposed(1,1)=Pmexp*4*rand(1,1); %probability of moving
        % thetaproposed(2,1)=Ppexp*4*rand(1,1); %probability of proliferating
        TSigma=2.38/sqrt(2)*[0.1,0.001];
        
        thetaProposed=thetaOld+TSigma'.*randn(2,1);
        thetaProposed(1)=(thetaProposed(1)>0.99)*abs(1.98-thetaProposed(1))+(thetaProposed(1)<0.99)*abs(thetaProposed(1));
        thetaProposed(2)=(thetaProposed(2)>0.01)*abs(0.02-thetaProposed(2))+(thetaProposed(2)<0.01)*abs(thetaProposed(2));
        
        % Check if we should accept or reject thetaProposed
        
        for reps=1:replicates
            %[Aproposed(:,:,reps),FRproposed(:,:,reps),FCproposed(:,:,reps)]=ABCmex(thetaProposed(1,1),thetaProposed(2,1),Mexp(:,:,reps));
            [Aproposed(:,:,reps),~,~]=ABCmex(thetaProposed(1,1),thetaProposed(2,1),Mexp(:,:,reps));
            %Bproposed(i)=(sum(sum(Aproposed(:,:,i))));
            % Compute SS: 1 #cells
            %S1proposed(i)=(sum(sum(Aproposed(:,:,i))));
            % Compute SS: 2 Trajectory distance
            %for n=1:N0
            %    Trajectory=[FRproposed(n,:,i);FCproposed(n,:,i)];
            %    F1(n,1)=sum(sum(abs(diff(Trajectory')')));
            %end
            %meanF1=mean(F1(1:N0,:),1); % average over the tracks, N0
            %SS6A=mean(meanF1); % average over the replicates
            %S2proposed(i)=SS6A;
            
            [I,J]=find(Aproposed(:,:,reps));
            D=mandist([I,J]'); 
            for rr=1:24
                Sproposedall(rr,reps)=length(find(triu(D)==rr));
            end
        end
        Sproposed=mean(Sproposedall,2);
       % Sproposed=[S1proposed;S2proposed];
        qhash=prod(exp(-diag((Sproposed-Sexp)'*Weights*(Sproposed-Sexp))./(24*2*Delta^2)));
        %qhash=prod(exp(-(Bproposed-Bexp).^2/(2*Delta^2)));
        alpha=min([qhash/qstar,1]);
        omega=rand(1);
        
        if omega>alpha
            theta(:,samples+1)=thetaOld;
            if p==P
                Sregression(:,samples+1)=Sregression(:,samples);
            end
        else
            theta(:,samples+1)=thetaProposed;
            thetaOld=thetaProposed;
            qstar=qhash;
            acceptedP(p)=acceptedP(p)+1;
            if p==P
                Sregression(:,samples+1)=Sproposed;
            end
        end 
    end
    thetaRecord=[thetaRecord,theta];
end
disp([num2str(toc),' Time to compute tolerance scheme'])


%% At this point apply regression adjustment on accepted samples to get a better covariance:
% I need a vector of summary statistics in the final population:
% S1(i);S2(i)



tic;
% %samples=500;numstats=1;numparams=1;
% %theta=rand(numparams,samples);
% %S=rand(numstats,samples)+5*theta;
% %Sregression;
% %Sexp=rand(numstats,1);
% %delta=2;
 numstats=size(Sregression,1);
% if numstats==1
%     S1r=Sregression(1,:);
%     Snormed=[S1r-Sexp(1)];
% elseif numstats==2
%     S1r=Sregression(1,:);
%     S2r=Sregression(2,:);
%     Snormed=[S1r-Sexp(1);S2r-Sexp(2)];
% elseif numstats==3
%     S1r=Sregression(1,:);
%     S2r=Sregression(2,:);
%     S3r=Sregression(3,:);
%     Snormed=[S1r-Sexp(1);S2r-Sexp(2);S3r-Sexp(3)];
% else
%     Snormed=[];
%     for numstatscounter=1:numstats
%         Sr=Sregression(numstatscounter,:);
%         Snormed=[Snormed;(Sr-Sexp(numstatscounter))];
%         %Snormed=[Snormed;(Sr-Sexp(numstatscounter))/Ssigma(numstatscounter)];
%     end
%     %error('Cannot handle more than 3 stats');
% end
% % % COMMENT: I've noticed that trying to adapt the weights to favour closer
% % % parameters renders beta (the regression gradient) close to zero, (due 
% % % to ignoring too many data points) so does little to narrow the parameter distribution
% %weights=3*(1-(sqrt(sum(Snormed.^2,1))./Delta).^2)/(4*Delta);
% Xdesignmatrix=[ones(deltaSamplesScheme(P),1),Snormed'];
% %Paramweights=ones(1,2000);
% Paramweights=ones(1,length(Sregression));
% %weights(weights<0)=0;
% 
% AB=(Xdesignmatrix'*diag(Paramweights)*Xdesignmatrix)^(-1) * Xdesignmatrix' *diag(Paramweights)* theta';
% Alpha=AB(1,:);
% Beta=AB(2:end,:);
% thetanew=(theta'-Snormed'*Beta)';

thetanew=RegressionAdjustmentFunction(theta,Sregression,Sexp,Ssigma,0);

theta=thetanew;
disp([num2str(toc),' Time to perform regression adjustment'])



tic;
for q=1:Q % perform ABC MCMC for every delta
   
    if q==1
         mu=mean(theta,2); TSigma=cov(theta'); R=chol(TSigma);
        thetaProposedVector= repmat(mu,1,clonesSamplesScheme(q)) + R*randn(2,clonesSamplesScheme(q));
        theta=zeros(2,clonesSamplesScheme(q));    
        %theta(:,1)=mean(theta,2);
        theta(:,1)=mu;        
        
    else
        TSigma=cov(theta'); R=chol(TSigma);
        thetaProposedVector= repmat(mu,1,clonesSamplesScheme(q)) + R*randn(2,clonesSamplesScheme(q));
        thetaTemp=theta(:,q-1);
        theta=zeros(2,clonesSamplesScheme(q));    
        theta(:,1)=thetaTemp;
        %Delta=deltaScheme(q);
        K=clonesScheme(q);
    end
    
    S=zeros(24,K);
    Sproposed=zeros(24,K);
    Sall=zeros(24,K,replicates);
    Sproposedall=zeros(24,K,replicates);
    % Generate Data from thetaStar and compute qstar
    for i=1:K
        for reps=1:replicates
            [A(:,:,i,reps),~,~]=ABCmex(thetaOld(1,1),thetaOld(2,1),Mexp(:,:,reps));
             [I,J]=find(A(:,:,i,reps));
             D=mandist([I,J]'); 
             for rr=1:24
                 Sall(rr,i,reps)=length(find(triu(D)==rr));
             end
        end
    end
    S=mean(Sall,3);

    SdiffMat=zeros(numstats,K);
    for kIter=1:K
        SdiffMat(:,kIter)=S(:,kIter)-Sexp;
    end
    
    qstar=prod(exp(-diag(SdiffMat'*Weights*SdiffMat)./(24*2*Delta^2)));
    
    
    for samples=1:(clonesSamplesScheme(q)-1)
        

        thetaProposed=thetaProposedVector(:,samples);
        %thetaProposed=thetaOld+Sigma'.*randn(2,1);
        thetaProposed(1)=(thetaProposed(1)>0.99)*abs(1.98-thetaProposed(1))+(thetaProposed(1)<0.99)*abs(thetaProposed(1));
        thetaProposed(2)=(thetaProposed(2)>0.01)*abs(0.02-thetaProposed(2))+(thetaProposed(2)<0.01)*abs(thetaProposed(2));
        
    for i=1:K
        for reps=1:replicates
            [Aproposed(:,:,i,reps),~,~]=ABCmex(thetaProposed(1,1),thetaProposed(2,1),Mexp(:,:,reps));
             [I,J]=find(Aproposed(:,:,i,reps));
             D=mandist([I,J]'); 
             for rr=1:24
                 Sproposedall(rr,i,reps)=length(find(triu(D)==rr));
             end
        end
    end
    Sproposed=mean(Sproposedall,3);

    % 
   % S=[S1;S2];
  %  Sexp=[S1exp;S2exp];
    %Weights=diag([sigma(1),sigma(2)]);
    SdiffMat=zeros(numstats,K);
    for kIter=1:K
        SdiffMat(:,kIter)=Sproposed(:,kIter)-Sexp;
    end
    
    %qstar=prod(exp(-diag((S-Sexp)'*Weights*(S-Sexp))./(2*Delta^2)));
    qhash=prod(exp(-diag(SdiffMat'*Weights*SdiffMat)./(24*2*Delta^2)));
    
 
        
    
%        Sproposed=[S1proposed;S2proposed];
%        qhash=prod(exp(-diag((Sproposed-Sexp)'*Weights*(Sproposed-Sexp))./(2*Delta^2)));
    
        alpha=min([qhash/qstar,1]);
        omega=rand(1);
        
        if omega>alpha
            theta(:,samples+1)=thetaOld;
        else
            theta(:,samples+1)=thetaProposed;
            thetaOld=thetaProposed;
            qstar=qhash;
            acceptedQ(q)=acceptedQ(q)+1;
        end 
    
    end
    thetaRecord=[thetaRecord,theta];
end
acceptedAll=[acceptedP,acceptedQ];
disp([num2str(toc),' Time to compute clones scheme'])
