clc;
clear;
addpath('./../distmesh/');
config;
%% Rectangle
[p,t]=distmesh2d(@fdRec,@huniform,scale,[-0,-0;1,1],...
    [0,0; 1,0; 0,1; 1,1]);

[C,NC,CMid,NCMid,NeighbourCells]=buildConnection(t);

phi=fdNotchedCircle(p);

save(['notchedCircle_' num2str(scale) '.mat'],'C','NC','CMid','NCMid','NeighbourCells','t','p','phi');

%% 2 circles in a block
phi=fd2Circles(p);

save(['2Circles_' num2str(scale) '.mat'],'C','NC','CMid','NCMid','NeighbourCells','t','p','phi');

%% U
clf;

[p,t]=distmesh2d(@fdU,@huniform,scale,[-0,-0;1,1],...
    [0,0; 1,0; 0,1; 1,1; 0.4,0.4; 0.6,0.4; 0.4,1; 0.6,1]);

[C,NC,CMid,NCMid,NeighbourCells]=buildConnection(t);

phi=fdUpperBlock(p);

save(['U_Region_' num2str(scale) '.mat'],'C','NC','CMid','NCMid','NeighbourCells','t','p','phi');

%% Functions

function d = fdRec( p )
d=drectangle(p, -0, 1, -0, 1);
end

function d = fdU( p )
dBlock=drectangle(p, -0, 1, -0, 1);
dUR=fdUpperBlock(p);
d=ddiff(dBlock, dUR);
end

function d = fdNotchedCircle( p )
gRec1=drectangle(p, 0.43, 0.57, 0.65, 0.99);
gCircle=dcircle(p, 0.5, 0.65, 0.25);
d=ddiff(gCircle, gRec1);
end

function d = fdUpperBlock( p )
d1=drectangle(p, 0.4, 0.6, 0.4, 1.1);
d2=drectangle(p, 0, 1, 1, 1.1);
d=dunion(d1,d2);
end

function d = fd2Circles( p )
d1=dcircle(p, 0.3, 0.3, 0.2);
d2=dcircle(p, 0.7, 0.7, 0.2);
d=dunion(d1,d2);
end