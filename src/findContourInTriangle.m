function [lineSeg1, lineSeg2] = findContourInTriangle(coordinates, phi)
% coordinates*: [x, y; x, y; x, y] coordinates of 3 vertices
% phi*: 3x1 value of phi on 3 vertices
% output: lineSeg*: [x, y] 2 ends of the line segment

if any(size(phi) ~= [3, 1])
    error('size(phi) shall be 3x1');
end


isEdge12Cutted= (phi(1)>=0 && phi(2) < 0) || (phi(1)<=0 && phi(2) > 0);
isEdge23Cutted= (phi(2)>=0 && phi(3) < 0) || (phi(2)<=0 && phi(3) > 0);
isEdge31Cutted= (phi(3)>=0 && phi(1) < 0) || (phi(3)<=0 && phi(1) > 0);


if (isEdge31Cutted && isEdge23Cutted)
    % check on Edge 13 and 23
    lineSeg1 = findZeroOnEdgeWrapper(1, 3, coordinates, phi);
    lineSeg2 = findZeroOnEdgeWrapper(2, 3, coordinates, phi);
elseif (isEdge12Cutted && isEdge23Cutted)
    % check on Edge 12 and 23
    lineSeg1 = findZeroOnEdgeWrapper(1, 2, coordinates, phi);
    lineSeg2 = findZeroOnEdgeWrapper(2, 3, coordinates, phi);
elseif (isEdge12Cutted && isEdge31Cutted)
    % check on Edge 12 and 13
    lineSeg1 = findZeroOnEdgeWrapper(1, 2, coordinates, phi);
    lineSeg2 = findZeroOnEdgeWrapper(1, 3, coordinates, phi);
else
    error('All phi have the same sign!')
end

end

function [R] = findZeroOnEdgeWrapper(t1, t2, coordinates, phi)
R = findZeroOnEdge(coordinates(t1,:),coordinates(t2,:),phi(t1),phi(t2));
end