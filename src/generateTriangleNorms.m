function [triangleOutNorms ,cornerNorm, areas] = generateTriangleNorms(p, t)
% input: p, t: mesh
% output: triangleOutNorm of size (nTriangles, 3, 2)
%         cornerNorm of size (nTriangles, 3, 2)

nTriangles=size(t, 1);
triangleOutNorms=zeros(nTriangles, 3, 2);
cornerNorm=zeros(nTriangles, 3, 2);
areas=zeros(nTriangles, 1);

for idxT=1:nTriangles
    triangleIdxP=t(idxT,:);
    triangleCoordinates=p(triangleIdxP,:);

    % edge out norm: e12, e23, e31
    TON=triangleOutNorm(triangleCoordinates);
    triangleOutNorms(idxT, :, :)=TON;
    
    % corner norm
    % Nc1= e12 + e31
    cornerNorm(idxT, 1, :)=normr(TON(1, :) + TON(3, :));
    % Nc2= e12 + e23
    cornerNorm(idxT, 2, :)=normr(TON(1, :) + TON(2, :));
    % Nc3= e23 + e31
    cornerNorm(idxT, 3, :)=normr(TON(2, :) + TON(3, :));
    
    areas(idxT, 1)=polyarea(triangleCoordinates(:,1), triangleCoordinates(:,2));
end
end

