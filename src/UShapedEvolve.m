clc;
clear;
config;
load(['U_Region_' num2str(scale) '.mat']);

NPoints=size(phi,1);
Rsf=0.5;

% set velocity
velocity=ones(size(p,1),1);
velocity(p(:,1)<0.6)=Rsf;
velocity(p(:,1)<0.15)=Rsf*Rsf;

disp('Building edges ...')
ovData=buildOutgoingEdges(p,C,NC,CMid,NCMid);
[A,edgeWeights]=buildMatrixA(ovData);
disp('Building edges done.')

disp('Building triangle basic data ...')
[triangleOutNorms ,cornerNorm, areas]=generateTriangleNorms(p, t);
disp('Building triangle basic data done.')

NEvolove=round(0.124/(evolveStep*Rsf*Rsf));
Result=cell(NEvolove,1);
AnalyticalResult=cell(NEvolove,1);
EmResult=zeros(NEvolove,1);
EaResult=zeros(NEvolove,1);

isoLines=[-0.2:0.02:-0.04, -0.02:0.01:0.02, 0.04:0.02:0.4];

for r=1:30
    disp(['Reinitial ' num2str(r) ' ' datestr(now,13)]);
    
    GReinitial=calcGradient(p, phi, ovData, A, phi, edgeWeights);
    %GReinitial=calcGradientViaIntegral(t, triangleOutNorms ,cornerNorm, areas, phi);
    S=phi./sqrt(phi.^2+(scale*0.8)^2);
    deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStep;
    phi=phi-deltaAmountReinitial;
end

e = 0;
for idx=1:NEvolove
    disp(['Step ' num2str(idx) ' ' datestr(now,13)]);
    
    %% Evolve
    GEvolve=calcGradient(p, phi, ovData, A, velocity, edgeWeights);
    GEvolveModule=sqrt(GEvolve(:,1).^2 + GEvolve(:,2).^2);
    deltaAmountEvolve=GEvolveModule .* velocity * evolveStep;
    intermediaPhi=phi-0.5*deltaAmountEvolve;
    
    GEvolve=calcGradient(p, intermediaPhi, ovData, A, velocity, edgeWeights);
    GEvolveModule=sqrt(GEvolve(:,1).^2 + GEvolve(:,2).^2);
    deltaAmountEvolve=GEvolveModule .* velocity * evolveStep;
    phi=phi-deltaAmountEvolve;
    
    %% Reinitial
    NR=3;
    for r=1:NR
        disp(['Reinitial ' num2str(r) ' ' datestr(now,13)]);
        
        GReinitial=calcGradient(p, phi, ovData, A, phi, edgeWeights);
        %GReinitial=calcGradientViaIntegral(t, triangleOutNorms ,cornerNorm, areas, phi);
        S=phi./sqrt(phi.^2+(scale*1)^2);
        deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStep;
        intermediaPhi=phi-0.5*deltaAmountReinitial;
        
        GReinitial=calcGradient(p, intermediaPhi, ovData, A, phi, edgeWeights);
        %GReinitial=calcGradientViaIntegral(t, triangleOutNorms ,cornerNorm, areas, intermediaPhi);
        S=phi./sqrt(intermediaPhi.^2+(scale*1)^2);
        deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStep;
        phi=phi-deltaAmountReinitial;
    end
    
    %% Error % record
    e = e + evolveStep*Rsf*Rsf;
    phiAnalytical = UShapedAnalytical(e, p, scale/100);
    
    Result{idx, 1} = phi;
    AnalyticalResult{idx, 1} = phiAnalytical;
    
    ErrorTags=abs(phiAnalytical)<5*scale;
    Error=abs(phi(ErrorTags)-phiAnalytical(ErrorTags));
    EmResult(idx)=max(Error);
    EaResult(idx)=mean(Error);
    disp(['Error: max->' num2str(EmResult(idx)) ' average->' num2str(EaResult(idx))]);
    
    hold on;
    [~, h]=tricontour(p,t,phi,[0,0]);
    for idxH=1:numel(h)
        h(idxH).EdgeColor='b';
    end
    [~, h]=tricontour(p,t,phiAnalytical,[0,0]);
    for idxH=1:numel(h)
        h(idxH).EdgeColor='r';
    end
    drawnow;
end

save(['resultData/UShaped_' num2str(scale) '.mat'], 'Result', 'AnalyticalResult', 'EmResult', 'EaResult', 'p', 't');
%save(['resultData/UShaped_ExplictDisabled_' num2str(scale) '.mat'], 'Result', 'AnalyticalResult', 'EmResult', 'EaResult', 'p', 't');