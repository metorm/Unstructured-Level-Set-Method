clc;
clear;
config;

%% load mesh
load(['3Circles_' num2str(scale) '.mat']);
NPoints=size(phi,1);

%% initial phi
x1=0.25;
x2=0.75;
y1=0.25;
y2=0.75;

phi=-min(min(min(-y1+p(:,2), y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));
phi=phi.*((p(:,1)-0.3).^2 + (p(:,2)-0.2).^2 + 0.1);

%% auxiliary data
disp('Building edges ...')
ovData=buildOutgoingEdges(p,C,NC,CMid,NCMid);
[A,edgeWeights]=buildMatrixA(ovData);
disp('Building edges done.')

%% evolve
NEvolove=80;
NWriteInterval=2;
Result=cell(round(NEvolove/NWriteInterval)+1,1);

isoLines=-0.5:1e-2:0.5;

for e=1:NEvolove
    disp(['Step ' num2str(e) ' ' datestr(now,13)]);
    
    if mod(e,NWriteInterval)==1
        Result{1+(e-1)/NWriteInterval,1}=phi;
    end
    
    GReinitial=calcGradient(p, phi, ovData, A, phi, edgeWeights);
    S=phi./sqrt(phi.^2+(scale*1)^2);
    deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStep;
    intermediaPhi=phi-0.5*deltaAmountReinitial;
    
    GReinitial=calcGradient(p, intermediaPhi, ovData, A, phi, edgeWeights);
    S=phi./sqrt(intermediaPhi.^2+(scale*1)^2);
    deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStep;
    phi=phi-deltaAmountReinitial;
    
end

%% draw
clf;
tricontour(p,t,phi,isoLines);
xlim([0,1]);
ylim([0,1]);
drawnow;