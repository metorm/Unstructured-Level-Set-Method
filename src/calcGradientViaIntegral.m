function [ G ] = calcGradientViaIntegral(t, triangleOutNorms ,cornerNorm, areas, phi)

nPoints=size(phi,1);
nTriangles=size(t, 1);

%% per triangle
trG=zeros(nTriangles, 2);

for idxT=1:nTriangles
    
    % out norm: e12, e23, e31
    TON=squeeze(triangleOutNorms(idxT, :, :));
    % phi at v1, v2, v3
    trianglePhi=phi(t(idxT,:));
    
    % gradient
    trG(idxT, :) = sum(TON .* (trianglePhi + trianglePhi([2, 3, 1]))) * 0.5 / areas(idxT);
end

%% per vertex
G=zeros(nPoints, 2);

for idxP=1:nPoints
    vertexTag=(t==idxP);
    neighbourCellTags=(any(vertexTag, 2));
    nNeighbourCells=sum(neighbourCellTags);

    neighbourCellGradient=trG(neighbourCellTags, :);
    crtCornerNorm=reshape(cornerNorm(neighbourCellTags, :, :), [3*nNeighbourCells, 2]);
    
    vertexTag=vertexTag(neighbourCellTags, :);
    vertexTag=reshape(vertexTag, [3*nNeighbourCells, 1]);
    crtCornerNorm(~vertexTag, :)=[];
    
    originWeight=dot(normr(neighbourCellGradient), crtCornerNorm, 2);
    weights=max(originWeight, 0);
    if sum(weights) < eps
        weights = originWeight - min(originWeight);
    end
    G(idxP, :) = sum(weights .* neighbourCellGradient, 1)/sum(weights);
end
end

