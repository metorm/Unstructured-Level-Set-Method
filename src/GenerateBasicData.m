clc;
clear;
addpath('./../distmesh/');
config;

[p,t]=distmesh2d(@fdRec,@huniform,scale,[-0,-0;1,1],...
    [-0,-0;1,-0;-0,1;1,1]);

[C,NC,CMid,NCMid]=buildConnection(t);

phi=fdNotchedCircle(p);

save(['notchedCircle_' num2str(scale) '.mat'],'C','NC','CMid','NCMid','t','p','phi');

%patch('vertices',p,'faces',t,'edgecol','k','FaceVertexCData',(phi>0)*1,'FaceColor','interp');

function d = fdRec( p )
d=drectangle(p, -0, 1, -0, 1);
end

function d = fdNotchedCircle( p )
gRec1=drectangle(p, 0.43, 0.57, 0.6, 0.99);
gCircle=dcircle(p, 0.5, 0.65, 0.25);
d=ddiff(gCircle, gRec1);
end