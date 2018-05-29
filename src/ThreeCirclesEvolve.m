clc;
clear;
config;

%% load mesh
load(['3Circles_' num2str(scale) '.mat']);
NPoints=size(phi,1);

%% initialization
circleR=0.2;
d1=dcircle(p, 0.3, 0.3, circleR);
d2=dcircle(p, 0.7, 0.3, circleR);
d3=dcircle(p, 0.5, 0.7, circleR);
phi=dunion(dunion(d1,d2),d3);

velocity=ones(size(p,1),1);

%% auxiliary data
disp('Building edges ...')
ovData=buildOutgoingEdges(p,C,NC,CMid,NCMid);
[A,edgeWeights]=buildMatrixA(ovData);
disp('Building edges done.')

%% evolve
NEvolove=100;
Result=cell(NEvolove,1);
EmResult=zeros(NEvolove,1);
EaResult=zeros(NEvolove,1);

isoLines=[-0.2:0.02:-0.04, -0.02:0.01:0.02, 0.04:0.02:0.4];

for e=1:NEvolove
    disp(['Step ' num2str(e) ' ' datestr(now,13)]);
    Result{e,1}=phi;
        
    GEvolve=calcGradient(p, phi, ovData, A, velocity, edgeWeights);
    GEvolveModule=sqrt(GEvolve(:,1).^2 + GEvolve(:,2).^2);
    deltaAmountEvolve=GEvolveModule .* velocity * evolveStep;
    intermediaPhi=phi-0.5*deltaAmountEvolve;
    
    GEvolve=calcGradient(p, intermediaPhi, ovData, A, velocity, edgeWeights);
    GEvolveModule=sqrt(GEvolve(:,1).^2 + GEvolve(:,2).^2);
    deltaAmountEvolve=GEvolveModule .* velocity * evolveStep;
    phi=phi-deltaAmountEvolve;
    
    %% error & record
    circleR=circleR+velocity * evolveStep;
    %% analytical data
    d1=dcircle(p, 0.3, 0.3, circleR);
    d2=dcircle(p, 0.7, 0.3, circleR);
    d3=dcircle(p, 0.5, 0.7, circleR);
    phiAnalytical=dunion(dunion(d1,d2),d3);
    
    ErrorTags=abs(phiAnalytical)<5*scale;
    Error=abs(phi(ErrorTags)-phiAnalytical(ErrorTags));
    EmResult(e)=max(Error);
    EaResult(e)=mean(Error);
    disp(['Error: max->' num2str(EmResult(e)) ' average->' num2str(EaResult(e))]);
end

%% save
save(['resultData/3Circles_' num2str(scale) '.mat'], 'Result', 'EmResult', 'EaResult', 'p', 't');