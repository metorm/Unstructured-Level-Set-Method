clc;
clear;
config;
load(['notchedCircle_' num2str(scale) '.mat']);

NPoints=size(phi,1);

centerBasedCoordinate=p-center;
velocity=normr(centerBasedCoordinate);
velocity=[velocity(:,2),-velocity(:,1)] .* sqrt(centerBasedCoordinate(:,1).^2+centerBasedCoordinate(:,2).^2);

disp('Building edges ...')
ovData=buildOutgoingEdges(p,C,NC,CMid,NCMid);
A=buildMatrixA(ovData);
disp('Building edges done.')


NEvolove=2000;
NWriteInterval=5;
Result=cell(round(NEvolove/NWriteInterval)+1,1);

for r=1:80
    disp(['Reinitial ' num2str(r) ' ' datestr(now,13)]);
    
    GReinitial=calcGradient(phi, ovData, A, phi);
    S=phi./sqrt(phi.^2+(scale*0.8)^2);
    deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStep;
    phi=phi-deltaAmountReinitial;
    
    if mod(r,5)==1
        clf;
        subplot(1,2,1);
        tricontour(p,t,phi,-0.2:0.05:0.2);
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
    
    GEvolve=calcGradient(phi, ovData, A, velocity);
    crtNormalVelocity=dot(GEvolve,velocity,2);
    %deltaAmountEvolve=crtNormalVelocity .* sqrt(GEvolve(:,1).^2+GEvolve(:,2).^2) * evolveStep;
    deltaAmountEvolve=crtNormalVelocity * evolveStep;
    phi=phi-deltaAmountEvolve;
    
    NR=6;
    for r=1:NR
        disp(['Reinitial ' num2str(r) ' ' datestr(now,13)]);
        
        GReinitial=calcGradient(phi, ovData, A, phi);
        S=phi./sqrt(phi.^2+(scale*0.8)^2);
        deltaAmountReinitial=S .* (sqrt(GReinitial(:,1).^2+GReinitial(:,2).^2)-1) * reinitialStep;
        phi=phi-deltaAmountReinitial;
    end
    
    if mod(e,5)==1
        clf;
        subplot(2,2,1);
        %trisurf(t,p(:,1),p(:,2),phi,'EraseMode','xor');
        %patch('vertices',p,'faces',t,'edgecol','k','FaceVertexCData',S,'FaceColor','interp','EraseMode','xor');
        tricontour(p,t,phi,-0.2:0.05:0.2);
        xlim([0,1]);
        ylim([0,1]);
        %view(40,130);
        subplot(2,2,2);
        %patch('vertices',p,'faces',t,'edgecol','k','FaceVertexCData',S,'FaceColor','interp','EraseMode','xor');
        quiver(p(:,1),p(:,2),GEvolve(:,1),GEvolve(:,2));
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
        %patch('vertices',p,'faces',t,'edgecol','k','FaceVertexCData',phi,'FaceColor','interp','EraseMode','xor');
        quiver(p(:,1),p(:,2),GReinitial(:,1),GReinitial(:,2))
        xlim([0,1]);
        ylim([0,1]);
        drawnow;
    end
end