function plot_control_vect(r,v,qi,wb,Ib,DT)

wbx = [ 0 -wb(3) wb(2);...
       wb(3) 0 -wb(1);...
       -wb(2) wb(1) 0];

%% get the forward, starboard, nadir quaternion

% forward is velocity vector
F = v/norm(v);

% nadir is -r
N = -r/norm(r);

% starboard is -FxN
S = -cross(F,N);

R_fsn = [F',S',N'];
q_fsn = Quaternion_fromRotationMatrix(R_fsn);

%% get the error between the current quaternion and the FSN axes

% want the error in the body frame
q_err = Quaternion_MULTIPLY(Quaternion_inv(qi),q_fsn);

% the axis of the quaternion defines the optimal w vector
if(q_err(1) < 0) q_err = -q_err; end
w_des_dir = q_err(2:4)';
q_err(1)

% scale the desired vector by how large the error is
max_w_des = 0.001;
w_des = w_des_dir * max_w_des * (1-abs(q_err(1)));   

control_torques = Ib*(w_des - wb')/DT - wbx * Ib * wb';

%% plot the vectors
quiver3(r(1),r(2),r(3),w_des_dir(1),w_des_dir(2),w_des_dir(3),...
    'autoscale','off','color','b',LineStyle='-');
quiver3(r(1),r(2),r(3),control_torques(1),control_torques(2),control_torques(3),...
    'autoscale','off','color','r',LineStyle='--');
