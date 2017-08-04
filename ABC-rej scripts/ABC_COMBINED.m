function ABC_COMBINED(Ratio)
%% Compute synthetic experiments: data and statistics

ROWS=24;
replicates=10;
DensityFraction=1;
%Ratio=100;
% cd ~/Documents/MATLAB/
% load('DifferentLatticeSpacings.mat', 'M24')
% load('ABCsims_Mnew_50replicates.mat', 'Mnew')
% load('ABCsims_Mnew2_50replicates.mat', 'Mnew2')
Mexp=zeros(24,32,replicates);
for i=1:replicates
%    %initialise in first 5 rows:
%    M=zeros(32,24);
%    M(randsample(32*5,DensityFraction*24))=1;
%    Mexp(:,:,i)=M';
     %initialise uniformly:
     M=zeros(24,32);
     M(randsample(32*24,DensityFraction*24))=1;
     Mexp(:,:,i)=M;
end

%% Step 0: Set up ABC parameters

%N0=sum(sum(Mexp));
N0=5;

if Ratio==10
    Pm=0.025;
    Pp=0.0025;
elseif Ratio==100
    Pm=0.25;
    Pp=0.0025;
elseif Ratio==1000
    Pm=0.25;
    Pp=0.00025;
end

Aexp=false(ROWS,ROWS*4/3,replicates); % Final matrices
%B=zeros(replicates,samples); % Final number of cells
%E1=zeros(24,replicates,samples); % Total distance moved by first 24 cells
FRexp=zeros(N0,36,replicates); % One trajectory with 36 points (the cell with index 0)
FCexp=zeros(N0,36,replicates); % One trajectory with 36 points (the cell with index 0)


for j=1:replicates
    [Aexp(:,:,j),FRexp(:,:,j),FCexp(:,:,j)]=ABCmex(Pm,Pp,Mexp(:,:,j));
end

Bexp=zeros(1,replicates);
for j=1:replicates
    Bexp(j)=(sum(sum(Aexp(:,:,j))));
end


%% Step 3b: compute other summary stats
Cexp=zeros(replicates,ROWS*7/3); %2D correlations
Gexp=zeros(2,replicates); %Gyration tensor eigenvalues
%E2sim=zeros(2,replicates,samples); % Distance moved according to trajectory (x+y, and sqrt(x^2+y^2))

F1exp=zeros(N0,replicates); % Distance moved by N0 cells according to trajectory (x+y)
F2exp=zeros(N0,replicates); % Distance moved by N0 cells according to trajectory (sqrt(x^2+y^2))
S1exp=zeros(N0,replicates); % Straightness Index according to trajectory (x+y)
S2exp=zeros(N0,replicates); % Straightness Index according to trajectory (sqrt(x^2+y^2))



for j=1:replicates
    [I,J]=find(Aexp(:,:,j));
    D=mandist([I,J]'); %manhattan distance, though could use boxdist.

    for i=1:ROWS*7/3
        Cexp(j,i)=length(find(triu(D)==i));
    end

    Stemp=zeros(2);
    N=Bexp(j);
    for i=1:length(I)
        for ii=1:length(I)
            Stemp(2,2)=Stemp(2,2)+(I(i)-I(ii))*(I(i)-I(ii))/(2*N^2);
            Stemp(1,2)=Stemp(1,2)+(I(i)-I(ii))*(J(i)-J(ii))/(2*N^2);
            Stemp(2,1)=Stemp(2,1)+(I(i)-I(ii))*(J(i)-J(ii))/(2*N^2);
            Stemp(1,1)=Stemp(1,1)+(J(i)-J(ii))*(J(i)-J(ii))/(2*N^2);
        end
    end
    Gexp(1:2,j)=eig(Stemp);

    for n=1:N0

        Trajectory=[FRexp(n,:,j);FCexp(n,:,j)];

        F1exp(n,j)=sum(sum(abs(diff(Trajectory')')));
        F2exp(n,j)=sum(sqrt(sum(diff(Trajectory')'.^2)));

        S1exp(n,j)=sum(sum(abs(diff(Trajectory')')))/sum(abs(Trajectory(:,end)-Trajectory(:,1)));
        S2exp(n,j)=sum(sqrt(sum(diff(Trajectory')'.^2)))/sqrt(sum((Trajectory(:,end)-Trajectory(:,1)).^2));
    end
    j
end
Cexp=Cexp';

%



replicates=10; %?
%Ratio=100;

if sum(sum(sum(Mexp==2)))~=0
    error('Invalid input matrices')
end

%sum(sum(sum(Mexp==2)))
%Mexp(find(Mexp==2))=1;

%% ABC pseudo-code

ROWS=24;
%% Step 0: Set up ABC parameters
samples=10000;    

%N0=sum(sum(Mexp));
N0=5;

A=false(ROWS,ROWS*4/3,replicates,samples); % Final matrices
%B=zeros(replicates,samples); % Final number of cells
%E1=zeros(24,replicates,samples); % Total distance moved by first 24 cells
FR=zeros(N0,36,replicates,samples); % One trajectory with 36 points (the cell with index 0)
FC=zeros(N0,36,replicates,samples); % One trajectory with 36 points (the cell with index 0)

theta=zeros(2,samples); %collect Pm and Pp

if Ratio==10
    Pmexp=0.025;
    Ppexp=0.0025;
elseif Ratio==100
    Pmexp=0.25;
    Ppexp=0.0025;
elseif Ratio==1000
    Pmexp=0.25;
    Ppexp=0.00025;
elseif Ratio==323
    Pmexp=0.025;
    Ppexp=0.00125;
end

%% Step 1: Draw parameter from priors
    %adapt priors depending on what you expect, for 231 cells, prior
    %Pm 0->0.05 and Pp 0->0.005
    % adapt priors so that Pm+Pp<=1
    if 4*(Pmexp+Ppexp)>1
        PmPpsum=4*(Pmexp+Ppexp);
        Pmexp=Pmexp/PmPpsum;
        Ppexp=Ppexp/PmPpsum;
    end
    theta(1,1:samples)=Pmexp*4*rand(1,samples); %probability of moving
    theta(2,1:samples)=Ppexp*4*rand(1,samples); %probability of proliferating

    
    
%% Step 2+3a: Simulate data from model with sampled parameters and collect summary stats S
tic;
parfor i=1:samples
    for j=1:replicates
        [A(:,:,j,i),FR(:,:,j,i),FC(:,:,j,i)]=ABCmex(theta(1,i),theta(2,i),Mexp(:,:,j));
    end
end
toc

B=zeros(replicates,samples);
for i=1:samples
    for j=1:replicates
        B(j,i)=(sum(sum(A(:,:,j,i))));
    end
end


%% Step 3b: compute other summary stats
Csim=zeros(replicates,ROWS*7/3,samples); %2D correlations
Gsim=zeros(2,replicates,samples); %Gyration tensor eigenvalues
%E2sim=zeros(2,replicates,samples); % Distance moved according to trajectory (x+y, and sqrt(x^2+y^2))

F1sim=zeros(N0,replicates,samples); % Distance moved by N0 cells according to trajectory (x+y)
F2sim=zeros(N0,replicates,samples); % Distance moved by N0 cells according to trajectory (sqrt(x^2+y^2))
S1sim=zeros(N0,replicates,samples); % Straightness Index according to trajectory (x+y)
S2sim=zeros(N0,replicates,samples); % Straightness Index according to trajectory (sqrt(x^2+y^2))



for j=1:replicates
    for k=1:samples
        [I,J]=find(A(:,:,j,k));
        D=mandist([I,J]'); %manhattan distance, though could use boxdist.
        
        for i=1:ROWS*7/3
            Csim(j,i,k)=length(find(triu(D)==i));
        end
        
        Stemp=zeros(2);
        N=B(j,k);
        for i=1:length(I)
            for ii=1:length(I)
                Stemp(2,2)=Stemp(2,2)+(I(i)-I(ii))*(I(i)-I(ii))/(2*N^2);
                Stemp(1,2)=Stemp(1,2)+(I(i)-I(ii))*(J(i)-J(ii))/(2*N^2);
                Stemp(2,1)=Stemp(2,1)+(I(i)-I(ii))*(J(i)-J(ii))/(2*N^2);
                Stemp(1,1)=Stemp(1,1)+(J(i)-J(ii))*(J(i)-J(ii))/(2*N^2);
            end
        end
        Gsim(1:2,j,k)=eig(Stemp);
        
        for n=1:N0

            Trajectory=[FR(n,:,j,k);FC(n,:,j,k)];

            F1sim(n,j,k)=sum(sum(abs(diff(Trajectory')')));
            F2sim(n,j,k)=sum(sqrt(sum(diff(Trajectory')'.^2)));

            S1sim(n,j,k)=sum(sum(abs(diff(Trajectory')')))/sum(abs(Trajectory(:,end)-Trajectory(:,1)));
            S2sim(n,j,k)=sum(sqrt(sum(diff(Trajectory')'.^2)))/sqrt(sum((Trajectory(:,end)-Trajectory(:,1)).^2));
        end
        
    end
    j
end

%save(['ABCmatsRatio_',int2str(Ratio),'_4_6_15.mat']);
%save('ABCmatsRatio_3T3_12_6_15.mat');

