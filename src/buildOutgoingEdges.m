function OV = buildOutgoingEdges(p, C, NC, CMid, NCMid)

NVertices=size(NC,1);
OV=cell(NVertices,4);

xPos=[1, 0];
xNeg=[-1, 0];
yPos=[0, 1];
yNeg=[0, -1];
divergenceThreshold=cos(deg2rad(50));

for idxP=1:NVertices
    crdT=[p(idxP, 1) p(idxP, 2)];
    
    neighborTag=C(idxP,1:NC(idxP))';
    crdV=p(neighborTag,:);
    crdOV=crdV-crdT;
    crdOVNorm=normr(crdOV);
    
    neighborTagMid=CMid(idxP,1:NCMid(idxP),:);
    neighborTagMid=reshape(neighborTagMid, [NCMid(idxP) 2]);
    crdVMid=(p(neighborTagMid(:,1),:)+p(neighborTagMid(:,2),:))./2;
    crdOVMid=crdVMid-crdT;
    crdOVNormMid=normr(crdOVMid);
    
    % xPos -> 1
    divergence=dot(crdOVNorm',repmat(xPos,NC(idxP),1)')';
    logicTag=divergence>divergenceThreshold;
    divergenceMid=dot(crdOVNormMid',repmat(xPos,NCMid(idxP),1)')';
    logicTagMid=divergenceMid>divergenceThreshold;
    OV{idxP,1}=compactOVInfo(logicTag, crdOV, neighborTag, logicTagMid, crdOVMid, neighborTagMid);
    
    % xNeg -> 2
    divergence=dot(crdOVNorm',repmat(xNeg,NC(idxP),1)')';
    logicTag=divergence>divergenceThreshold;
    divergenceMid=dot(crdOVNormMid',repmat(xNeg,NCMid(idxP),1)')';
    logicTagMid=divergenceMid>divergenceThreshold;
    OV{idxP,2}=compactOVInfo(logicTag, crdOV, neighborTag, logicTagMid, crdOVMid, neighborTagMid);
    
    % yPos -> 3
    divergence=dot(crdOVNorm',repmat(yPos,NC(idxP),1)')';
    logicTag=divergence>divergenceThreshold;
    divergenceMid=dot(crdOVNormMid',repmat(yPos,NCMid(idxP),1)')';
    logicTagMid=divergenceMid>divergenceThreshold;
    OV{idxP,3}=compactOVInfo(logicTag, crdOV, neighborTag, logicTagMid, crdOVMid, neighborTagMid);
    
    % yNeg -> 4
    divergence=dot(crdOVNorm',repmat(yNeg,NC(idxP),1)')';
    logicTag=divergence>divergenceThreshold;
    divergenceMid=dot(crdOVNormMid',repmat(yNeg,NCMid(idxP),1)')';
    logicTagMid=divergenceMid>divergenceThreshold;
    OV{idxP,4}=compactOVInfo(logicTag, crdOV, neighborTag, logicTagMid, crdOVMid, neighborTagMid);
    
end

end

function OVInfoCell = compactOVInfo(logicTag, crdOV, neighborTag, logicTagMid, crdOVMid, neighborTagMid)
selectedcrdOV=crdOV(logicTag,:);
selectedneighborTag=neighborTag(logicTag);
selectedcrdOVMid=crdOVMid(logicTagMid,:);
selectedneighborTagMid=neighborTagMid(logicTagMid,:);
OVInfoCell={[selectedcrdOV, sqrt(selectedcrdOV(:,1).^2+selectedcrdOV(:,2).^2)], selectedneighborTag, ...
    [selectedcrdOVMid, sqrt(selectedcrdOVMid(:,1).^2+selectedcrdOVMid(:,2).^2)], selectedneighborTagMid};
end
