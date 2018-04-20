clc;
clear;
config;
load(['2Circles_' num2str(scale) '.mat']);

NPoints=size(phi,1);
Rsf=0.5;

velocity=ones(size(p,1),1);
velocity(p(:,1) - 0.5 < 0.5 - p(:,2))=Rsf;

disp('Building edges ...')
ovData=buildOutgoingEdges(p,C,NC,CMid,NCMid);
[A,edgeWeights]=buildMatrixA(ovData);
disp('Building edges done.')

NEvolove=80;
NWriteInterval=2;
Result=cell(round(NEvolove/NWriteInterval)+1,1);

isoLines=[-0.2:0.02:-0.04, -0.02:0.01:0.02, 0.04:0.02:0.4];

for r=1:30
    disp(['Reinitial ' num2str(r) ' ' datestr(now,13)]);
    
    GReinitial=calcGradient(p, phi, ovData, A, phi, edgeWeights, true);
    S=phi./sqrt(phi.^2+(scale*0.8)^2);
    deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStep;
    phi=phi-deltaAmountReinitial;
    
    if mod(r,1)==0
        clf;
        subplot(1,2,1);
        tricontour(p,t,phi,isoLines);
        xlim([0,1]);
        ylim([0,1]);
        subplot(1,2,2);
        quiver(p(:,1),p(:,2),GReinitial(:,1),GReinitial(:,2));
        xlim([0,1]);
        ylim([0,1]);
        drawnow;
    end
end

for e=1:NEvolove
    disp(['Step ' num2str(e) ' ' datestr(now,13)]);
    
    if mod(e,NWriteInterval)==1
        Result{1+(e-1)/NWriteInterval,1}=phi;
    end
    
    
    GEvolve=calcGradient(p, phi, ovData, A, velocity, edgeWeights, true);
    GEvolveModule=sqrt(GEvolve(:,1).^2 + GEvolve(:,2).^2);
    deltaAmountEvolve=GEvolveModule .* velocity * evolveStep;
    intermediaPhi=phi-0.5*deltaAmountEvolve;
    
    GEvolve=calcGradient(p, intermediaPhi, ovData, A, velocity, edgeWeights, true);
    GEvolveModule=sqrt(GEvolve(:,1).^2 + GEvolve(:,2).^2);
    deltaAmountEvolve=GEvolveModule .* velocity * evolveStep;
    phi=phi-deltaAmountEvolve;
    
    NR=10;
    for r=1:NR
        disp(['Reinitial ' num2str(r) ' ' datestr(now,13)]);
        
        GReinitial=calcGradient(p, phi, ovData, A, phi, edgeWeights, true);
        S=phi./sqrt(phi.^2+(scale*1)^2);
        deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStep;
        intermediaPhi=phi-0.5*deltaAmountReinitial;
        
        GReinitial=calcGradient(p, intermediaPhi, ovData, A, phi, edgeWeights, true);
        S=phi./sqrt(intermediaPhi.^2+(scale*1)^2);
        deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStep;
        phi=phi-deltaAmountReinitial;
        
    end
    
    if mod(e,1)==0
        clf;
        subplot(2,2,1);
        grid on;
        %trisurf(t,p(:,1),p(:,2),phi,'EraseMode','xor');
        %patch('vertices',p,'faces',t,'edgecol','k','FaceVertexCData',S,'FaceColor','interp','EraseMode','xor');
        tricontour(p,t,phi,isoLines);
        xlim([0,1]);
        ylim([0,1]);
        %view(40,130);
        subplot(2,2,2);
        quiver(p(:,1),p(:,2),GEvolve(:,1),GEvolve(:,2));
        %patch('vertices',p,'faces',t,'LineStyle','none',...
        %    'FaceVertexCData',sqrt(GEvolve(:,1).^2+GEvolve(:,2).^2),...
        %    'FaceColor','interp');
        xlim([0,1]);
        ylim([0,1]);
        subplot(2,2,3);
        %trisurf(t,p(:,1),p(:,2),sqrt(G(:,1).^2+G(:,2).^2),'EraseMode','xor');
        %trisurf(t,p(:,1),p(:,2),deltaAmountEvolve,'EraseMode','xor');
        tricontour(p,t,deltaAmountEvolve,10);
        xlim([0,1]);
        ylim([0,1]);
        %view(40,130);
        subplot(2,2,4);
        %patch('vertices',p,'faces',t,'LineStyle','none',...
        %    'FaceVertexCData',sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2),...
        %    'FaceColor','interp');
        quiver(p(:,1),p(:,2),GReinitial(:,1),GReinitial(:,2))
        xlim([0,1]);
        ylim([0,1]);
        drawnow;
    end
end

save resultData/2Circles_0.025_Explicit.mat Result