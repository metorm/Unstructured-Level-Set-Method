res=100;
center=[0.5, 0.5];

x=linspace(0,1,res);
y=x;
[Gx, Gy]=meshgrid(x,y);
u=(pi/3.14)*(center(2)-Gy);
v=(pi/3.14)*(Gx-center(1));

p=[reshape(Gx, [res*res, 1]), reshape(Gy, [res*res, 1])];
centerBasedCoordinate=p-center;
velocity=normr(centerBasedCoordinate);
velocity=[-velocity(:,2),velocity(:,1)] .* sqrt(centerBasedCoordinate(:,1).^2+centerBasedCoordinate(:,2).^2);

hold on;
quiver(Gx,Gy,u,v);
quiver(p(:,1), p(:,2), velocity(:,1), velocity(:,2));
