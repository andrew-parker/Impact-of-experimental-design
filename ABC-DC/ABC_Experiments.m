%% Make experimental data independent of the ABC_DC or ABC_rej algorithms:
replicates=10;
DensityFraction=DF;
ScratchDensity=SD;
%% 1: Compute initial domains

Pm=0.25;
Pp=0.0025;
ROWS=24;
Mexp=zeros(24,32,replicates);
for reps=1:replicates
M=zeros(32,24);
M(randsample(32*ScratchDensity,DensityFraction*24))=1;
Mexp(:,:,reps)=M';
end
%load(['sigma_SD',num2str(SD),'DF',num2str(DF),'.mat'])
%Ssigma=sigmalarge(2:25);

Aexp=false(ROWS,ROWS*4/3,replicates); % Final matrices cropped

%% 2: Simulate to find experimental SS
Sexpall=zeros(24,replicates);
for reps=1:replicates
    [Aexp(:,:,reps),~,~]=ABCmex(Pm,Pp,Mexp(:,:,reps));
    [I,J]=find(Aexp(:,:,reps));
    D=mandist([I,J]'); %manhattan distance, though could use boxdist.
    for rr=1:24
        Sexpall(rr,reps)=length(find(triu(D)==rr));
    end
end
Sexp=mean(Sexpall,2);
%%
