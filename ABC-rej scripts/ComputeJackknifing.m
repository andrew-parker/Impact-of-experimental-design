%% To get an approximate stddev for IKL, want to perform some resampling on the IKLs, which involved recomputing distances, and then recomputing the information gain from the combination we are interested in!
%% First: Compute distance 10 times (once for each replicate) and store Dlargeresampled (:,:,1:10)
% Initialise, and compute statistics before they are averaged over the
% replicates:
replicates=10;
N0=5;
samples=20000;
DlargeRE=zeros(20,samples,replicates);
tic
SS6Aexpinit=mean(F1exp(1:N0,:),1); %size reps * 10, average over the 10 tracks for size = reps,1
SS6Asiminit=mean(F1sim(1:N0,:,:),1); %size=N0, reps, samples mean over first 10 tracks
% Trajectory-calculated average distance moved by 5 cells (dx^2+dy^2)^0.5
SS6Bexpinit=mean(F2exp(1:N0,:),1); %size reps * 10, average over the 10 tracks for size = reps,1
SS6Bsiminit=mean(F2sim(1:N0,:,:),1); %size=N0, reps, samples mean over first 10 tracks

SS7expinitA=zeros(1,replicates);
SS7expinitB=zeros(1,replicates);
SS7initA=zeros(replicates,samples);
SS7initB=zeros(replicates,samples);
for j=1:replicates
    for i=1:samples
        SS7initA(j,i)=mean(S1sim(isinf(S1sim(:,j,i))==0,j,i));
        SS7initB(j,i)=mean(S2sim(isinf(S2sim(:,j,i))==0,j,i));
    end
    SS7expinitA(j)=mean(S1exp(isinf(S1exp(:,j))==0,j));
    SS7expinitB(j)=mean(S2exp(isinf(S2exp(:,j))==0,j));
end



%% Compute new statistics
L1=zeros(replicates,samples);
L1exp=zeros(1,replicates);
L2=zeros(replicates,samples);
L2exp=zeros(1,replicates);
%% Clustersize:
disp(['Computing clustersize'])
for j=1:replicates
    %disp(['Computed clustersize for replicate ', int2str(j)])
    %EXP
    CC=bwconncomp(Aexp(:,:,j),8);
    L1exp(j)=max(cellfun('length',CC.PixelIdxList));
    CC=bwconncomp(Aexp(:,:,j),4);
    L2exp(j)=max(cellfun('length',CC.PixelIdxList));
    %SIM
    
    for i=1:samples
        CC=bwconncomp(A(:,:,j,i),8);
        L1(j,i)=max(cellfun('length',CC.PixelIdxList));
        CC=bwconncomp(A(:,:,j,i),4);
        L2(j,i)=max(cellfun('length',CC.PixelIdxList));
    end
end

%% QuadratCounts:

H=zeros(replicates,samples,3);
Hexp=zeros(replicates,3);

for widthIndex=1:3
    widthofBoxes=2^widthIndex; % must be 8, or 4 or 2 (some other common divisor of 24 and 32)
    num_bX=24/widthofBoxes;
    num_bY=32/widthofBoxes;
    num_boxes=num_bX*num_bY;

    disp(['Computing quadrats for widthIndex',int2str(widthIndex)])
for j=1:replicates
    %disp(['Computed quadrats for replicate ', int2str(j)])
    %EXP
    n=zeros(1,num_boxes);
    nbar=Bexp(j)/num_boxes;
    for ii=1:num_bX
        for jj=1:num_bY
            n(jj+num_bY*(ii-1))=sum(sum(Aexp(widthofBoxes*(ii-1)+1:widthofBoxes*ii,widthofBoxes*(jj-1)+1:widthofBoxes*(jj),j)));
        end
    end
    Hexp(j,widthIndex)=sum((n-nbar).^2);
    
    %SIM
    for i=1:samples
        n=zeros(1,num_boxes);
        nbar=B(j,i)/num_boxes;
        for ii=1:num_bX
            for jj=1:num_bY
                n(jj+num_bY*(ii-1))=sum(sum(A(widthofBoxes*(ii-1)+1:widthofBoxes*ii,widthofBoxes*(jj-1)+1:widthofBoxes*(jj),j,i)));
            end
        end
        H(j,i,widthIndex)=sum((n-nbar).^2);
    end
end
end

%% Pair Correlations
%compute Ctilde
ROWS=24;
COLS=32;
Mfull=ones(ROWS,COLS);
[I,J]=find(Mfull);
Dmat=mandist([I,J]');
Ctilde=zeros(1,ROWS+COLS);
for i=1:(ROWS+COLS)
Ctilde(i)=length(find(triu(Dmat)==i));
end

%compute new correlations (assume we already have C)
C2exp=zeros(size(Cexp));
C2sim=zeros(size(Csim));
for j=1:replicates
    C2exp(:,j)=(Cexp(:,j)./(Ctilde'.*Bexp(j).*(Bexp(j)-1)./(ROWS*COLS*(ROWS*COLS-1))));
    for i=1:samples
        C2sim(j,:,i)=(Csim(j,:,i)./(Ctilde.*B(j,i).*(B(j,i)-1)./(ROWS*COLS*(ROWS*COLS-1))));
    end
end

%1D correlations
DDxexp=zeros(ROWS,replicates);
DDyexp=zeros(COLS,replicates);
Qxexp=zeros(ROWS-1,replicates);
Qyexp=zeros(COLS-1,replicates);

for j=1:replicates
    [I,J]=find(Aexp(:,:,j));
    N=Bexp(j);
    for k=1:N
        for m=(k+1):N
            DDxexp(abs(I(k)-I(m))+1,j)=DDxexp(abs(I(k)-I(m))+1,j)+1;
            DDyexp(abs(J(k)-J(m))+1,j)=DDyexp(abs(J(k)-J(m))+1,j)+1;
        end
    end
    Qxexp(:,j)=DDxexp(2:24,j)./(32*32*(24.-(1:23)')*(N/(24*32))*((N-1)/(24*32-1)));
    Qyexp(:,j)=DDyexp(2:32,j)./(24*24*(32.-(1:31)')*(N/(24*32))*((N-1)/(24*32-1)));
end

DDxsim=zeros(ROWS,replicates,samples);
DDysim=zeros(COLS,replicates,samples);
Qxsim=zeros(ROWS-1,replicates,samples);
Qysim=zeros(COLS-1,replicates,samples);

disp(['Computing pairwise corr'])
for j=1:replicates
    %disp(['Computed pairwise corr for replicate ',int2str(j)])
    for i=1:samples
        [I,J]=find(A(:,:,j,i));
        N=B(j,i);
        for k=1:N
            for m=(k+1):N
                DDxsim(abs(I(k)-I(m))+1,j,i)=DDxsim(abs(I(k)-I(m))+1,j,i)+1;
                DDysim(abs(J(k)-J(m))+1,j,i)=DDysim(abs(J(k)-J(m))+1,j,i)+1;
            end
        end
        Qxsim(:,j,i)=DDxsim(2:24,j,i)./(32*32*(24.-(1:23)')*(N/(24*32))*((N-1)/(24*32-1)));
        Qysim(:,j,i)=DDysim(2:32,j,i)./(24*24*(32.-(1:31)')*(N/(24*32))*((N-1)/(24*32-1)));
    end
end
toc





for replicateExcluded=1:10
    X=[1:replicateExcluded-1,replicateExcluded+1:10];
    disp('Computing old statistics')
    tic
           
    D2c=zeros(1,samples);
    D4=zeros(1,samples);
%% Compute Old statistics

%% Compute trajectory averaged statistics

SS1exp=mean(Bexp(X),2);
SS1sim=mean(B(X,:),1); %length of this vector = #samples

SS2exp=mean(Cexp(:,X),2); % correlations normalised size 10,56 -> 1,56
SS2sim=zeros(ROWS*7/3,samples);
for i=1:samples
    SS2sim(:,i)=mean(Csim(X,:,i),1); % Csim has dim 10,56,100, mean(Csim,1) has dim 1,56,100. 
end
%DEPENDS ON Cexp ORIENTATION, SHOULD BE (mean(CMexp(:,X),2))
SS3exp=max(mean(Cexp(X,:),1)); % Peak height of non-normalised data
SS3sim=squeeze(max(mean(Csim(X,:,:),1)));
SS4exp=mean(Cexp(:,X),2)/sum(mean(Cexp(:,X),2)); % correlations normalised size 10,56 -> 1,56
SS4sim=zeros(ROWS*7/3,samples);
for i=1:samples
    SS4sim(:,i)=mean(Csim(X,:,i),1)/sum(mean(Csim(X,:,i),1)); % Csim has dim 10,56,100, mean(Csim,1) has dim 1,56,100. 
end
%SS4exp=mean(E1exp(1,X),2); % Total distance moved by one cell
%SS4exp=mean(mean(Eexp(1:5,:))); %Total distance moved by 5 cells
%SS4sim=squeeze(mean(E1(1,X,:),2));
%SS4sim=squeeze(mean(mean(E(1:5,:,:))));
SS5Aexp=mean(Gexp(1,X),2);
SS5Asim=squeeze(mean(Gsim(1,X,:),2));
SS5Bexp=mean(Gexp(2,X),2);
SS5Bsim=squeeze(mean(Gsim(2,X,:),2));
%SS5exp=SS5Aexp+SS5Bexp; % Sum of Eigenvalues of Gyration Tensor
%SS5sim=SS5Asim+SS5Bsim;
SS6Aexp=mean(SS6Aexpinit(:,X),2); % Trajectory-calculated distance moved by 1 cell (dx+dy)
SS6Asim=squeeze(mean(SS6Asiminit(:,X,:),2));
SS6Bexp=mean(SS6Bexpinit(:,X),2); % Trajectory-calculated distance moved by 1 cell (dx+dy)
SS6Bsim=squeeze(mean(SS6Bsiminit(:,X,:),2));
SS7Aexp=mean(SS7expinitA(:,X),2);
SS7Bexp=mean(SS7expinitB(:,X),2);
SS7Asim=mean(SS7initA(X,:),1);
SS7Bsim=mean(SS7initB(X,:),1);



D1c=abs(SS1sim-SS1exp)/mad(SS1sim,1);
sigma2c=mad(SS2sim(1:24,:)',1)';
sigma4=mad(SS4sim(1:24,:)',1)';
for i=1:samples
    D2c(i)=(sum(((SS2sim(1:24,i)-SS2exp(1:24))./sigma2c).^2)/24).^.5;
    %D4(i)=sum(abs(SS4sim(:,i)-SS4exp));
    D4(i)=(sum(((SS4sim(1:24,i)-SS4exp(1:24))./sigma4).^2)/24).^.5;
end
D3=abs(SS3sim-SS3exp)'/mad(SS3sim,1);
D5A=abs(SS5Asim-SS5Aexp)'/mad(SS5Asim,1);
D5B=abs(SS5Bsim-SS5Bexp)'/mad(SS5Bsim,1);
D6A=abs(SS6Asim-SS6Aexp)'/mad(SS6Asim,1);
D6B=abs(SS6Bsim-SS6Bexp)'/mad(SS6Bsim,1);
D7A=abs(SS7Asim-SS7Aexp)/mad(SS7Asim,1);
D7B=abs(SS7Bsim-SS7Bexp)/mad(SS7Bsim,1);

disp('Computing new statistics')

%% Computing new statistics distances


% Average distance
% D8A=zeros(1,samples);
% D8B=zeros(1,samples);
D8A=abs(mean(L1(X,:),1)-mean(L1exp(1,X),2))/mad(mean(L1(X,:),1),1);
D8B=abs(mean(L2(X,:),1)-mean(L2exp(1,X),2))/mad(mean(L2(X,:),1),1);

% D9A=zeros(1,samples);
% D9B=zeros(1,samples);
% D9C=zeros(1,samples);
D9A=abs(mean(H(X,:,1),1)-mean(Hexp(X,1),1))/mad(mean(H(X,:,1),1),1);
D9B=abs(mean(H(X,:,2),1)-mean(Hexp(X,2),1))/mad(mean(H(X,:,2),1),1);
D9C=abs(mean(H(X,:,3),1)-mean(Hexp(X,3),1))/mad(mean(H(X,:,3),1),1);


SSsim10A=zeros(24,samples);
SSexp10A=mean(C2exp(1:24,X),2); %2d real normed
SSsim10B=zeros(24,samples);
SSexp10B=mean(DDxexp(:,X),2); %1d X
SSsim10C=zeros(32,samples);
SSexp10C=mean(DDyexp(:,:),2); %1d Y
SSsim10D=zeros(23,samples);
SSexp10D=mean(Qxexp(:,:),2); %1d X normed
SSsim10E=zeros(31,samples);
SSexp10E=mean(Qyexp(:,:),2); %1d Y normed


for i=1:samples
   % SSsim(:,i)=mean(Csim(:,:,i),1); % Csim has dim 10,56,100, mean(Csim,1) has dim 1,56,100. 
   % SSsim(:,i)=mean(Csim(:,:,i),1)/sum(mean(Csim(:,:,i),1)); % Csim has dim 10,56,100, mean(Csim,1) has dim 1,56,100. 
    SSsim10A(:,i)=mean(C2sim(X,1:24,i),1); %2d real normed
    SSsim10B(:,i)=mean(DDxsim(:,X,i),2);%1d X
    SSsim10C(:,i)=mean(DDysim(:,X,i),2);%1d Y
    SSsim10D(:,i)=mean(Qxsim(:,X,i),2);%1d X normed
    SSsim10E(:,i)=mean(Qysim(:,X,i),2);%1d Y normed
end

sigma10A=mad(SSsim10A(1:24,:)',1)';
sigma10B=mad(SSsim10B(1:24,:)',1)';
sigma10C=mad(SSsim10C(1:32,:)',1)';
sigma10D=mad(SSsim10D(1:23,:)',1)';
sigma10E=mad(SSsim10E(1:31,:)',1)';

D10A=zeros(1,samples);
D10B=zeros(1,samples);
D10C=zeros(1,samples);
D10D=zeros(1,samples);
D10E=zeros(1,samples);

for i=1:samples
    D10A(i)=(sum(((SSsim10A(1:24,i)-SSexp10A(1:24))./sigma10A).^2)/24).^.5;
    D10B(i)=(sum(((SSsim10B(1:24,i)-SSexp10B(1:24))./sigma10B).^2)/24).^.5;
    D10C(i)=(sum(((SSsim10C(1:32,i)-SSexp10C(1:32))./sigma10C).^2)/32).^.5;
    D10D(i)=(sum(((SSsim10D(1:23,i)-SSexp10D(1:23))./sigma10D).^2)/23).^.5;
    D10E(i)=(sum(((SSsim10E(1:31,i)-SSexp10E(1:31))./sigma10E).^2)/31).^.5;
end

    
DlargeRE(:,:,replicateExcluded)=[D1c;D2c;D3;D4;D5A;D5B;D6A;D6B;D7A;D7B;D8A;D8B;D9A;D9B;D9C;D10A;D10B;D10C;D10D;D10E]; %20 rows x 20000 cols x 10 replicates

toc
end


%% Now we have DlargeRE, we compute information for each replicateRemoved:

% First run: just compute for all statistics...
STC=1:20; % StatsToConsider
COMBOS1=nchoosek(STC,1);
COMBOS2=nchoosek(STC,2);
COMBOS3=nchoosek(STC,3);

accCI1=[1,11:15,7,9,5,2,16,17,19]'; %no 3,4,6,8,10,18,20
accCI2=find(ones(length(COMBOS2),1)-(COMBOS2(:,1)==3 | COMBOS2(:,1)==4 | COMBOS2(:,1)==6 | COMBOS2(:,1)==8 | COMBOS2(:,1)==10 | COMBOS2(:,1)==18 | COMBOS2(:,1)==20 | COMBOS2(:,2)==3 | COMBOS2(:,2)==4 | COMBOS2(:,2)==6 | COMBOS2(:,2)==8 | COMBOS2(:,2)==10 | COMBOS2(:,2)==18 | COMBOS2(:,2)==20));
accCI3=find(ones(length(COMBOS3),1)-(COMBOS3(:,1)==3 | COMBOS3(:,1)==4 | COMBOS3(:,1)==6 | COMBOS3(:,1)==8 | COMBOS3(:,1)==10 | COMBOS3(:,1)==18 | COMBOS3(:,1)==20 | COMBOS3(:,2)==3 | COMBOS3(:,2)==4 | COMBOS3(:,2)==6 | COMBOS3(:,2)==8 | COMBOS3(:,2)==10 | COMBOS3(:,2)==18 | COMBOS3(:,2)==20 | COMBOS3(:,3)==3 | COMBOS3(:,3)==4 | COMBOS3(:,3)==6 | COMBOS3(:,3)==8 | COMBOS3(:,3)==10 | COMBOS3(:,3)==18 | COMBOS3(:,3)==20 ));

information1RE=zeros(length(COMBOS1),3,10);
information2RE=zeros(length(COMBOS2),3,10);
information3RE=zeros(length(COMBOS3),3,10);

for replicateExcluded=1:10

Dlarge=DlargeRE(:,:,replicateExcluded);
    
information1=zeros(length(COMBOS1),3);
information2=zeros(length(COMBOS2),3);
information3=zeros(length(COMBOS3),3);

for HowManyCombos=1:3
    disp(['Computing for HowManyCombos=',int2str(HowManyCombos)])
    switch HowManyCombos
        case 1
            COMBOS=COMBOS1;
            accCI=accCI1;
            information=zeros(size(COMBOS,1),3);
        case 2
            COMBOS=COMBOS2;
            accCI=accCI2;
            information=zeros(size(COMBOS,1),3);
        case 3
            COMBOS=COMBOS3;
            accCI=accCI3;
            information=zeros(size(COMBOS,1),3);
    end
%for kk=1:size(COMBOS,1)
for kk=1:size(accCI,1)
%for kk=1:1 %test
STC=COMBOS(accCI(kk),:); %statistics to consider are the rows of COMBOS#

D=sum(Dlarge(STC,:).^2,1);
[~,I]=sort(D);

%Compute posterior:
k=200;
%Epsilon(i,1)=D(I(k));
data=[theta(1,I(1:k))',theta(2,I(1:k))'];
n=512;
[bandwidth,density,X,Y]=kde2d(data,n,[0,0],[1,0.01]); % Note: relective density at boundary (0,1) and (0, 0.01)
%contourf(X,Y,density,'LineStyle','None');
density=density/trapz(X(1,:),trapz(Y(:,1),density,1));
prior=ones(n)*100;
density(density<=0)=1e-25;
density1=trapz(Y(:,1),density,1);
density2=trapz(X(1,:),density,2);
density1(density1<=0)=1e-25;
density2(density2<=0)=1e-25;
prior1=trapz(Y(:,1),prior,1);
prior2=trapz(X(1,:),prior,2);
Info2D=trapz(Y(:,1),trapz(X(1,:),density.*log(abs(density)./prior)));
Info1DA=trapz(X(1,:),(density1.*log(abs(density1)./prior1)));
Info1DB=trapz(Y(:,1),(density2.*log(abs(density2)./prior2)));
pause(0.01)
information(accCI(kk),1:3)=[Info2D,Info1DA,Info1DB];    
    
end
    switch HowManyCombos
        case 1
            information1=information;
        case 2
            information2=information;
        case 3
            information3=information;
    end
end
information1RE(:,:,replicateExcluded)=information1;
information2RE(:,:,replicateExcluded)=information2;
information3RE(:,:,replicateExcluded)=information3;
end

save(['infoRE_SD',int2str(ScratchDensity),'DF',int2str(DensityFraction),'.mat'],'DlargeRE','information1RE','information2RE','information3RE','theta');


