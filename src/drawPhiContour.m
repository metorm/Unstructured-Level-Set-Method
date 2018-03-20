clf;

conValues=[-0.2, -1e-2, 0, 1e-2, 0.2];

for idxP=1:size(Result)
    crtPhi=Result{idxP};
    if size(crtPhi,1) < 1 || size(crtPhi,1) < 2
        break;
    end
    
    pause(0.1);
    
    clf;
    tricontour(p,t,crtPhi,conValues);
    xlim([0,1]);
    ylim([0,1]);
    drawnow;
end