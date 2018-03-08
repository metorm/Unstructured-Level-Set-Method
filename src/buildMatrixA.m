function [rA] = buildMatrixA(ovData)
NPoints=size(ovData,1);
rA(NPoints,4).OK=true;
rA(NPoints,4).A=[];

for idxP=1:NPoints
    for idxC=1:4
        OV=ovData{idxP,idxC}{1,1};
        OVMid=ovData{idxP,idxC}{1,3};
        A=[OV(:,1:2);OVMid(:,1:2)];
        rA(idxP,idxC).A=A;
        rA(idxP,idxC).OK=size(A,1)>=2;
    end
end
end
