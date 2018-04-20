function [lineSeg1, lineSeg2] = findContourInTriangle(coordinates, phi)
% coordinates*: [x, y; x, y; x, y] coordinates of 3 vertices
% phi*: 3x1 value of phi on 3 vertices
% output: lineSeg*: [x, y] 2 ends of the line segment

if sum(size(phi) == [3, 1]) ~= 2
    error('size(phi) shall be 3x1');
end

signArray = sign(phi);

if signArray(1) == signArray(2) && signArray(1) ~= signArray(3)
    % check on Edge 13 and 23
    lineSeg1 = findZeroOnEdgeWrapper(1, 3, coordinates, phi);
    lineSeg2 = findZeroOnEdgeWrapper(2, 3, coordinates, phi);
elseif signArray(1) == signArray(3) && signArray(1) ~= signArray(2)
    % check on Edge 12 and 23
    lineSeg1 = findZeroOnEdgeWrapper(1, 2, coordinates, phi);
    lineSeg2 = findZeroOnEdgeWrapper(2, 3, coordinates, phi);
elseif signArray(2) == signArray(3) && signArray(2) ~= signArray(1)
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