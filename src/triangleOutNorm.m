function [TNorms] = triangleOutNorm(triangleCoordinates)
% input: [x1, y1; x2, y2; x3; y3]
% output: [n_12_x, n_12_y; n_23_x, n_23_y; n_31_x, n_31_y]

assert(all(size(triangleCoordinates)==[3, 2]));

triangleCoordinates=triangleCoordinates-mean(triangleCoordinates);
triangleCoordinatesShift=triangleCoordinates([2,3,1],:);
edgeCenters=(triangleCoordinates+triangleCoordinatesShift)./2;

edges=triangleCoordinates-triangleCoordinatesShift;
TNorms=edges(:,[2,1]).*[1 -1];

tagsToFlip=dot(edgeCenters, TNorms, 2) < 0;
TNorms(tagsToFlip, :) = TNorms(tagsToFlip, :) * -1;
end

