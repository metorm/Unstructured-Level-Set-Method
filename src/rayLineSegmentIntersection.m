function [Intersetion, L, isCrossed] = rayLineSegmentIntersection(...
    rayOrigin, rayDirection, lineSegEnd1, lineSegEnd2)
% input: [x,y] rayOrigin, rayDirection, lineSegEnd1, lineSegEnd2
% output: Intersetion: [x,y] corss point
%         L: length from rayOrigin to Intersetion
%         isCrossed: true/false indicating if corss point is found

rayDirection=normr(rayDirection);
lineSegVector=lineSegEnd2-lineSegEnd1;
segDirection=normr(lineSegVector);
segLength=norm(lineSegVector);

% if they are parallel
if abs(dot(rayDirection,segDirection)-1)<eps
    isCrossed=false;
    Intersetion=NaN;
    L=NaN;
    return;
end

% parametric equation
% rayOrigin + rt * rayDirection = lineSegEnd1 + st * segDirection

rtst=[rayDirection',-segDirection']\((lineSegEnd1-rayOrigin)');
rt=rtst(1);
st=rtst(2);

isCrossed= rt>0 && st>0 && st<segLength;
if isCrossed
    Intersetion=rayOrigin+rt*rayDirection;
    L=rt;
else
    Intersetion=NaN;
    L=NaN;
end

end