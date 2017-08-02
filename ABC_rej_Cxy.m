function [theta,I2c]=ABC_rej_Cxy(SD,DF)
ROWS=24;
replicates=10;
load(['ABCexperimentalCXYSD',int2str(SD),'DF',int2str(DF),'.mat'])
load(['sigma_SD',num2str(SD),'DF',num2str(DF),'.mat'])
Ssigma=sigmalarge(2:25);

SS2exp=Sexp;
sigma2c=Ssigma;
%% Step 0: Compute synthetic experiments: data and statistics
% % Ratio=100;
% % ROWS=24;
% % replicates=10;
% % %DensityFraction=1;
% % %Ratio=100;[
% % % cd ~/Documents/MATLAB/
% % % load('DifferentLatticeSpacings.mat', 'M24')
% % % load('ABCsims_Mnew_50replicates.mat', 'Mnew')
% % % load('ABCsims_Mnew2_50replicates.mat', 'Mnew2')
% % Mexp=zeros(24,32,replicates);
% % for i=1:replicates
% %    %initialise in first ScratchDensity rows:
% %    M=zeros(32,24);
% %    M(randsample(32*ScratchDensity,DensityFraction*24))=1;
% %    Mexp(:,:,i)=M';
% % %      %initialise uniformly:
% % %      M=zeros(24,32);
% % %      M(randsample(32*24,DensityFraction*24))=1;
% % %      Mexp(:,:,i)=M;
% % end
% % 
% % %% Step 0a: Set up ABC parameters
% % 
% % 
% % if Ratio==10
% %     Pm=0.025;
% %     Pp=0.0025;
% % elseif Ratio==100
% %     Pm=0.25;
% %     Pp=0.0025;
% % elseif Ratio==1000
% %     Pm=0.25;
% %     Pp=0.00025;
% % end
% % 
% % Aexp=false(ROWS,ROWS*4/3,replicates); % Final matrices
% % 
% % tic;
% % for j=1:replicates
% % %     [Aexp(:,:,j),FRexp(:,:,j),FCexp(:,:,j)]=ABCmex(Pm,Pp,Mexp(:,:,j));
% %     [Aexp(:,:,j),~,~]=ABCmex(Pm,Pp,Mexp(:,:,j));
% % end
% % disp([num2str(toc),' Time to compute experiments'])
% 
% %% Step 0b: compute experimental summary stats
% Cexp=zeros(replicates,ROWS*7/3); %2D correlations
% tic;
% for j=1:replicates
%     [I,J]=find(Aexp(:,:,j));
%     D=mandist([I,J]'); %manhattan distance, though could use boxdist.
%     for i=1:ROWS*7/3
%         Cexp(j,i)=length(find(triu(D)==i));
%     end
% end
% Cexp=Cexp';
% disp([num2str(toc),' Time to compute experimental statistics'])
% 
% %% Step 1: ABC code (after computing experimental statistics)
% 
% 
% if sum(sum(sum(Mexp==2)))~=0
%     error('Invalid input matrices')
% end

%% Step 1a: Set up ABC parameters

samples=20000;     %Changed to 20000

A=false(ROWS,ROWS*4/3,replicates,samples); % Final matrices

theta=zeros(2,samples); %collect Pm and Pp

%% Step 1b: Draw parameters from priors
    %adapt priors depending on what you expect, for 231 cells, prior
    %Pm 0->0.05 and Pp 0->0.005
    % adapt priors so that Pm+Pp<=1
    PmMax=0.99;
    PpMax=0.01;
    theta(1,1:samples)=PmMax*rand(1,samples); %probability of moving
    theta(2,1:samples)=PpMax*rand(1,samples); %probability of proliferating

    
    
%% Step 2: Simulate data from model with sampled parameters
tic;
for i=1:samples
    for j=1:replicates
       % [A(:,:,j,i),FR(:,:,j,i),FC(:,:,j,i)]=ABCmex(theta(1,i),theta(2,i),Mexp(:,:,j));
        [A(:,:,j,i),~,~]=ABCmex(theta(1,i),theta(2,i),Mexp(:,:,j));
    end
end
disp([num2str(toc),' Time to do simulations'])

%% Steo 3a: Collect summary stats S from A, FR, FC

%% Step 3b: compute other summary stats
Csim=zeros(replicates,ROWS*7/3,samples); %2D correlations
tic;
for j=1:replicates
    for k=1:samples
        [I,J]=find(A(:,:,j,k));
        D=mandist([I,J]'); %manhattan distance, though could use boxdist.
        
        for i=1:ROWS*7/3
            Csim(j,i,k)=length(find(triu(D)==i));
        end
    end
    %disp(['Calculated statistics for replicate ',int2str(j)]);
end
disp([num2str(toc),' Time to compute simulation statistics'])

X=1:replicates;
D2c=zeros(1,samples);
%SS2exp=mean(Cexp(:,X),2); % correlations normalised size 10,56 -> 1,56

SS2sim=zeros(ROWS*7/3,samples);
for i=1:samples
    SS2sim(:,i)=mean(Csim(X,:,i),1); % Csim has dim 10,56,100, mean(Csim,1) has dim 1,56,100. 
end

%sigma2c=mad(SS2sim(1:24,:)',1)';
for i=1:samples
    D2c(i)=(sum(((SS2sim(1:24,i)-SS2exp(1:24))./sigma2c).^2)/24).^.5;
end

[~,I2c]= sort(D2c); %corr
