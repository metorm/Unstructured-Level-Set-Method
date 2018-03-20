function [rA, edgeWeights] = buildMatrixA(ovData)
NPoints=size(ovData,1);
rA(NPoints,4).OK=true;
rA(NPoints,4).A=[];
edgeWeights=cell(NPoints,4);

for idxP=1:NPoints
    for idxC=1:4
        OV=ovData{idxP,idxC}{1,1};
        OVMid=ovData{idxP,idxC}{1,3};
        A=[OV(:,1:2);OVMid(:,1:2)];
        
        edgeWeights{idxP,idxC}=[ovData{idxP,idxC}{1,5};ovData{idxP,idxC}{1,6}];
        
        rA(idxP,idxC).A=A .* edgeWeights{idxP,idxC};
        rA(idxP,idxC).OK=size(A,1)>=2;
    end
end
end
