function  plot_space_craft(quat,trans)


% axis lim scales 
axlim = 1;
spacecraft_color = [0 0.5 0.5];
% rotation matrix 
dcm = RotationMatrix_from_quaternion(quat);

% generate a surface representing a spacecraft -- cylinder
fig = figure(1); clf(fig); hold on;
set(fig,'units','inches','position',[0 0 13 11],'color','w');

% specs
r = 0.03;
l = 0.2;

len_pts     = linspace(-l,l,10);
thet_pts    = linspace(0,2*pi,10);
r_pts       = linspace(0,r,3);

% create tube body
[L,T]   = meshgrid(len_pts,thet_pts);
X       = r.*cos(T);
Y       = r.*sin(T);
Z       = L;
[X,Y,Z] = rotateData(X,Y,Z,dcm);
[X,Y,Z] = translateData(X,Y,Z,trans);
tube = surf(X,Y,Z,'FaceColor',spacecraft_color);

% create lid
[R,T]   = meshgrid(r_pts,thet_pts);
X       = R.*cos(T);
Y       = R.*sin(T);
Z       = l*ones(size(X,1),size(X,2));
[X,Y,Z] = rotateData(X,Y,Z,dcm);
[X,Y,Z] = translateData(X,Y,Z,trans);
lid = surf(X,Y,Z,'FaceColor',spacecraft_color);

% create bottom
[R,T]   = meshgrid(r_pts,thet_pts);
X       = R.*cos(T);
Y       = R.*sin(T);
Z       = -l*ones(size(X,1),size(X,2));
[X,Y,Z] = rotateData(X,Y,Z,dcm);
[X,Y,Z] = translateData(X,Y,Z,trans);
bottom = surf(X,Y,Z,'FaceColor',spacecraft_color);

% plot direction arrows
f =  4*r*dcm*[1;0;0];
r =  4*r*dcm*[0;1;0];
d =  3/2*l*dcm*[0;0;1];
 
quiver3(trans(1),trans(2),trans(3),f(1),f(2),f(3),'autoscale','off',color='r');
quiver3(trans(1),trans(2),trans(3),r(1),r(2),r(3),'autoscale','off','color',[0 0.5 0]);
quiver3(trans(1),trans(2),trans(3),d(1),d(2),d(3),'autoscale','off',color='b');

axis equal; box on; xlim(axlim*[-1 1]); ylim(axlim*[-1 1]); zlim(axlim*[-1 1]);
view([30,20])

end

%% functions
function [X,Y,Z] = rotateData(X,Y,Z,R)

Rot_data = (R*[X(:),Y(:),Z(:)]')';
X = reshape(Rot_data(:,1),size(X,1),size(X,2));
Y = reshape(Rot_data(:,2),size(X,1),size(X,2));
Z = reshape(Rot_data(:,3),size(X,1),size(X,2));

end


function [X,Y,Z] = translateData(X,Y,Z,trans)

X = X + trans(1);
Y = Y + trans(2);
Z = Z + trans(3);

end