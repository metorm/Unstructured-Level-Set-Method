% this file contains necessary functions to calculate the U-shaped example
function [phi]=UShapedAnalytical(e, p, arcStep)
% main function
% input: e: propagate distance of the slowest area
%        p: points

%% coordinates of key vertices, left -> right
Rsf=0.5;
Rsfs2 = Rsf^2;
taper = pi/2 + asin(Rsf);
v1X = 0;
v1Y = 1 - e;

v3X = 0.15;
v3Y = 1 - e/Rsf;
[l23a, l23b, l23c] = pointK2ABC(v3X, v3Y, taper);

% line l12: y = v1Y
[v2X, v2Y] = lineCross(l23a, l23b, l23c, 0, 1, -v1Y);

edges = [...
    v1X, v1Y;...
    v2X, v2Y;...
    v3X, v3Y;...
    ];

% if v2X < 0, v1 disappears
if v2X < 0
    edges(1, :) = [];
    edges(1, :) = [0, -l23c/l23b;];
end

if (0.4 - e/Rsf) < 0.15
    error('Analytical solution fails when (0.4 - e/Rsf) < 0.15');
end
v3aX = 0.4 - e/Rsf;
v3aY = v3Y;

v4X = v3aX;
v4Y = 0.6;

edges = [edges; v3aX, v3aY; v4X, v4Y];

% there is a arc v4 ~ v5

v6X = 0.6;
v6Y = 0.6 - e/Rsfs2;
[l56a, l56b, l56c] = pointK2ABC(v6X, v6Y, taper);

v5aY = 0.6 - e/Rsf;
v5aX = (-l56c - l56b*v5aY) / l56a;

if v5aX > 0.4
    % horizon line v5 - v5a exists
    edges = [edges; arc2Poly(0.4, 0.6, e/Rsf, pi, 1.5*pi, arcStep);...
        v5aX, v5aY; v6X, v6Y];
else
    % the line starts from v6 intersects with arc v4 ~ v5
    [v5X, v5Y] = lineCrossCircle(l56a, l56b, l56c, 0.4, 0.6, e/Rsf, v6X, v6Y);
    theta2 = atan2(v5X - 0.4, v5Y - 0.6);
    while theta2 < pi
        theta2 = theta2 + 2*pi;
    end
    edges = [edges; arc2Poly(0.4, 0.6, e/Rsf, pi, theta2, arcStep); v6X, v6Y];
end

% there is a arc v6 ~ v7

v7X = 0.6 + e/Rsfs2;
v7Y = 0.6;

v9X = 1;
v9Y = 1 - e/Rsfs2;

v8X = v7X;
v8Y = v9Y;

if v8X > 1
    theta2 = 2*pi - acos(0.4/(e/Rsfs2));
    edges = [edges; arc2Poly(0.6, 0.6, e/Rsfs2, 1.5*pi, theta2, arcStep)];
else
   edges = [edges; arc2Poly(0.6, 0.6, e/Rsfs2, 1.5*pi, 2*pi, arcStep)]; 
   edges = [edges; v7X, v7Y; v8X, v8Y; v9X, v9Y];
end

%% distance
[~, phi, ~] = distance2curve(edges, p);

%% sign

if v2X > 0
    tag = p(:, 1) < v2X & p(:,2) > 1-e;
    phi(tag) = -phi(tag);
    tag = p(:, 1) >= v2X & p(:, 1) < v3X & isOnUpperOrRightOfLine(l23a, l23b, l23c, p(:,1), p(:,2));
    phi(tag) = -phi(tag);
else
    tag = p(:, 1) < v3X & isOnUpperOrRightOfLine(l23a, l23b, l23c, p(:,1), p(:,2));
    phi(tag) = -phi(tag);
end

tag = p(:, 1) >= v3X & p(:, 1) < v3aX & p(:, 2) > 1-e/Rsf;
phi(tag) = -phi(tag);

if v5aX > 0.4
    % horizon line v5 - v5a exists
    tag = p(:, 1) >= v3aX & p(:, 1) < 0.4 & ...
        (p(:, 2) >= 0.6 | ...
        (p(:, 2) < 0.6 & sqrt( (p(:, 1) - 0.4).^2 + (p(:, 2) - 0.6).^2) < e/Rsf));
    tag = tag | (p(:, 1) >= 0.4 & p(:, 1) < v5aX & p(:, 2) > 0.6 - e/Rsf);
    tag = tag | (p(:, 1) >= v5aX & p(:, 1) < 0.6 & ...
        isOnUpperOrRightOfLine(l56a, l56b, l56c, p(:,1), p(:,2)));
    
    phi(tag) = -phi(tag);
else
    % the line starts from v6 intersects with arc v4 ~ v5
    tag = p(:, 1) >= v3aX & p(:, 1) < v5X & ...
        (p(:, 2) >= 0.6 | ...
        (p(:, 2) < 0.6 & sqrt( (p(:, 1) - 0.4).^2 + (p(:, 2) - 0.6).^2) < e/Rsf));
    tag = tag | (p(:, 1) >= v5X & p(:, 1) < 0.6 & ...
        isOnUpperOrRightOfLine(l56a, l56b, l56c, p(:,1), p(:,2)));
    
    phi(tag) = -phi(tag);
end

if v8X > 1
    arcEndY = 0.6 - sqrt((e/Rsfs2)^2 - 0.16);
    tag = p(:, 1) >= 0.6 & (p(:, 2) >= arcEndY | ...
        (p(:, 2) < arcEndY & sqrt( (p(:, 1) - 0.6).^2 + (p(:, 2) - 0.6).^2) < e/Rsfs2));
    phi(tag) = -phi(tag);
else
    tag = p(:, 1) >= 0.6 & p(:, 1) < v7X & ...
        (p(:, 2) >= 0.6 | ...
        (p(:, 2) < 0.6 & sqrt( (p(:, 1) - 0.6).^2 + (p(:, 2) - 0.6).^2) < e/Rsfs2));
    tag = tag | (p(:, 1) >= v7X & p(:, 2) > 1-e/Rsfs2);
    
    phi(tag) = -phi(tag);
end

end

%% utility functions
function [A, B, C] = pointK2ABC(x0, y0, theta)
% converte line equation format
% (x0, y0), theta ==>> Ax + By + C = 0
%
% y-y0=tan(theta)(x-x0)
% y/tan(theta) - y0/tan(theta) = x - x0

iPi_2=theta/(pi/2);
iPi=theta/(pi);
isPi_2= abs(iPi_2-round(iPi_2)) < eps && abs(iPi - round(iPi)) >= eps;
isPi = abs(iPi - round(iPi))<eps;

if isPi_2
    A = 1;
    B = 0;
elseif isPi
    A = 0;
    B = 1;
else
    A = 1;
    B = -1/tan(theta);
end
C = - B*y0 - A*x0;
end

function [r] = isOnUpperOrRightOfLine(A, B, C, x, y)
if abs(B) < eps
    r = x > (-C/A);
else
    r = y > (-C - A*x)/B;
end
end

function [rX, rY] = lineCross(A1, B1, C1, A2, B2, C2)
r=[A1, B1; A2, B2]\[-C1;-C2];
rX=r(1);
rY=r(2);
end

function [rX, rY] = lineCrossCircle(A, B, C, cX, cY, R, xN, yN)
% find the intersection of circle and line, then choose the nearer one to
% [xN, yN] and return it
%
% Ax + By + C = 0
% x^2 - 2 cX  x + y^2 - 2 cY y + cX^2 + cY^2 = R^2
% (x - cX)^2 + (y - cY)^2 = R^2
if abs(A)<eps
    % y == -C/B
    yy = -C/B;
    t = sqrt(R^2-(yy - cY)^2);
    
    rX1 = cX-t;
    rY1 = yy;
    rX2 = cX+t;
    rY2 = yy;
elseif abs(B)<eps
    % x == -C/A
    xx = -C/A;
    t=sqrt(R^2-(xx - cX)^2);
    
    rX1 = xx;
    rY1 = cY-t;
    rX2 = xx;
    rY2 = cY+t;
else
    % y == -C/B - A/B x
    % y - cY = (-C/B - cY) - A/B x => t1 x + t2
    t1 = -A/B;
    t2 = -C/B - cY;
    % x^2 - 2 cX  x + cX^2 + t1^2 x^2 + 2 t1 t2 x + t2^2 = R^2
    % (1+t1^2) x^2 + (2 t1 t2 - 2cX) x + cX^2 + t2^2 - R^2 = 0
    b = 2 * t1 * t2 - 2 * cX;
    a = 1+t1^2;
    c = cX^2 + t2^2 - R^2;
    
    dlt = sqrt(b^2 - 4 * a * c);
    if ~isreal(dlt)
        if abs(imag(dlt)) < eps
            dlt = real(dlt);
        else
            error('No real root!');
        end
    end
    rX1 = (-b + dlt) / (2*a);
    rY1 = - (C + A*rX1) / B;
    rX2 = (-b - dlt) / (2*a);
    rY2 = - (C + A*rX2) / B;
end

d1 = (rX1 - xN)^2 + (rY1 - yN)^2;
d2 = (rX2 - xN)^2 + (rY2 - yN)^2;

if d1 > d2
    rX = rX2;
    rY = rY2;
else
    rX = rX1;
    rY = rY1;
end

end

function [r] = arc2Poly(cX, cY, R, theta1, theta2, stepScale)
% theta1, theta2 must be counter-clockwise

step = stepScale / R; % (2*pi) / (2*R*pi/stepScale)
thetas = (theta1:step:theta2)';
if abs(thetas(end) - theta2) > eps
    thetas = [thetas; theta2];
end

r = repmat([cX, cY], numel(thetas), 1) + [cos(thetas), sin(thetas)] * R;

end