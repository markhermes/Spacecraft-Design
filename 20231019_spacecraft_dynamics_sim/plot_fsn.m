function plot_fsn(r,v)

% forward is velocity vector
F = v/norm(v);

% nadir is -r
N = -r/norm(r);

% starboard is -FxN
S = -cross(F,N);

R_fsn = [F',S',N'];
q_fsn = Quaternion_fromRotationMatrix(R_fsn);

f = 0.5*Quaternion_rotateVectorWithQuat([0 1 0 0],q_fsn);
s = 0.5*Quaternion_rotateVectorWithQuat([0 0 1 0],q_fsn);
n = 0.5*Quaternion_rotateVectorWithQuat([0 0 0 1],q_fsn);

quiver3(r(1),r(2),r(3),f(1),f(2),f(3),'autoscale','off','color',[0.9,0.9,0.9],LineStyle='--');
quiver3(r(1),r(2),r(3),s(1),s(2),s(3),'autoscale','off','color',[0.9,0.9,0.9],LineStyle='--');
quiver3(r(1),r(2),r(3),n(1),n(2),n(3),'autoscale','off','color',[0.9,0.9,0.9],LineStyle='--');
end