clc;
clear;
config;
load(['RateStick_' num2str(scaleStick) '.mat']);

%% load mesh
NPoints=size(phi,1);

%% initialization
circleR=scaleStick/2;
phi=dcircle(p, 0.0, 0.025, circleR);

velocity=ones(size(p,1),1);

%% auxiliary data
disp('Building edges ...')
ovData=buildOutgoingEdges(p,C,NC,CMid,NCMid);
[A,edgeWeights]=buildMatrixA(ovData);
disp('Building edges done.')

%% evolve
NEvolove=round(1/evolveStepStick);
Result=cell(NEvolove,1);
EmResult=zeros(NEvolove,1);
EaResult=zeros(NEvolove,1);

isoLines=[-0.2:0.02:-0.04, -0.02:0.005:0.02, 0.04:0.02:0.4];

for e=1:NEvolove
    disp(['Step ' num2str(e) ' ' datestr(now,13)]);
    Result{e,1}=phi;
        
    GEvolve=calcGradient(p, phi, ovData, A, velocity, edgeWeights);
    GEvolveModule=sqrt(GEvolve(:,1).^2 + GEvolve(:,2).^2);
    deltaAmountEvolve=GEvolveModule .* velocity * evolveStepStick;
    intermediaPhi=phi-0.5*deltaAmountEvolve;
    
    GEvolve=calcGradient(p, intermediaPhi, ovData, A, velocity, edgeWeights);
    GEvolveModule=sqrt(GEvolve(:,1).^2 + GEvolve(:,2).^2);
    deltaAmountEvolve=GEvolveModule .* velocity * evolveStepStick;
    phi=phi-deltaAmountEvolve;
    
    %% error & record
    circleR=circleR + evolveStepStick;
    %% analytical data
    phiAnalytical=dcircle(p, 0.0, 0.025, circleR);
    
    ErrorTags=abs(phiAnalytical)<5*scaleStick;
    Error=abs(phi(ErrorTags)-phiAnalytical(ErrorTags));
    EmResult(e)=max(Error);
    EaResult(e)=mean(Error);
    disp(['Error: max->' num2str(EmResult(e)) ' average->' num2str(EaResult(e))]);
    
    NR=10;
    for r=1:NR
        disp(['Reinitial ' num2str(r) ' ' datestr(now,13)]);
        
        GReinitial=calcGradient(p, phi, ovData, A, phi, edgeWeights);
        S=phi./sqrt(phi.^2+(scaleStick*1)^2);
        deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStepStick;
        intermediaPhi=phi-0.5*deltaAmountReinitial;
        
        GReinitial=calcGradient(p, intermediaPhi, ovData, A, phi, edgeWeights);
        S=phi./sqrt(intermediaPhi.^2+(scaleStick*1)^2);
        deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStepStick;
        phi=phi-deltaAmountReinitial;
    end
end

%% save
save(['resultData/RateStick_' num2str(scaleStick) '.mat'], 'Result', 'EmResult', 'EaResult', 'p', 't');