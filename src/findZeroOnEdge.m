function [Z] = findZeroOnEdge(a,b,phiA,phiB)
% input: a&b: [x,y] line segment phiA/B: phi value on a and b
% output: Z: [x,y] point on ab where phi=0, or NaN if there's no zero point

if sign(phiA) == sign(phiB)
    Z=NaN;
    return;
end

Z=a+(abs(phiA)/(abs(phiA)+abs(phiB)))*(b-a);
end

