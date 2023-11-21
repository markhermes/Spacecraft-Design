function plot_earth()

r = 0.2;

thet    = linspace(-80,80,10)*pi/180;
phi     = linspace(0, 2*pi, 20);

[T,P]   = meshgrid(thet,phi);

Z = r*sin(T);
X = r*cos(T).*cos(P);
Y = r*cos(T).*sin(P);

surf(X,Y,Z,'FaceColor','b');

quiver3(0,0,0,2*r,0,0,'autoscale','off',color='r');
quiver3(0,0,0,0,2*r,0,'autoscale','off',color='g');
quiver3(0,0,0,0,0,2*r,'autoscale','off',color='b');


% also plot the equatorial plane
X = 5*r*cos(P);
Y = 5*r*sin(P);
Z = 0*X;

fill(X,Y,Z,'FaceColor','k','FaceAlpha',0.2);

end