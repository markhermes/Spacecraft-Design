% regulate the spacecraft so that the axes align with the Forward,
% Starboard, Nadir (FSN) axes

function [dx] = orbital_sim_with_fine_control(t,x,ge,Ib,Ib_inv,DT)

%% controller parameters
control_torques_max = 0.05; 

%% define variables
r = x(1:3);
v = x(4:6);
qb = x(7:10); 
wb = x(11:13);
wbx = [ 0 -wb(3) wb(2);...
       wb(3) 0 -wb(1);...
       -wb(2) wb(1) 0];

%% first translation -- simple

dr = v;
dv = -ge * r * (r'*r)^-(3/2);

%% get the forward, starboard, nadir quaternion

% forward is velocity vector
F = v/norm(v);

% nadir is -r
N = -r/norm(r);

% starboard is -FxN
S = -cross(F,N);

R_fsn = [F,S,N];
q_fsn = Quaternion_fromRotationMatrix(R_fsn);

%% apply the fine control

% want the error in the body frame
q_err = Quaternion_MULTIPLY(q_fsn,Quaternion_inv(qb));
q_err = q_err/norm(q_err);


% % integrate the angular value for the torque correction... Direction should
% % always be applied in the error direction (assuming small angular
% % velocities)
% angular_error_mag = acos(abs(q_err(1)))/2 ;
% angular_error_dir = (q_err(2:4)'/norm(q_err(2:4)) );
% de = angular_error_mag; 
% e_val = angular_error_mag;
% e_int = x(14);
% disp(sum(angular_error_mag));

%yaw pitch roll = z,y,x error
euler_err   = EulerAngle_fromQuaternionData_NED(q_err);
e_val       = flip(euler_err)';
e_int       = x(14:16);
de          = e_val;

% disp(norm(e_val));

% maybe first just apply torques in the w_des direction?
Kp = 1e0;
Ki = 5e-3;
Kd = 0.001;
pi_torque = (Kp*e_val + Ki*e_int);

% do the normal control in addition to the other torque
if(q_err(1) < 0) q_err = -q_err; end
w_des_dir = q_err(2:4)';
q_err(1);

% scale the desired vector by how large the error is
max_w_des = 0.2;
w_des = w_des_dir * max_w_des * (1.5-abs(q_err(1)));

% maybe first just apply torques in the w_des direction?
control_torques = pi_torque + (Ib*(w_des - wb)/DT - wbx * Ib * wb);
    

% set max value of control_torques
if(any(abs(control_torques) > control_torques_max))
    control_torques = control_torques_max* ...
        control_torques/max(abs(control_torques));
end

%% Then body rotation
% cannot find a nice solution for inertial rotation because the inertial
% tensor is not constant and would be a pain in the ass to calculate for
% each rotation

M = control_torques;

dwb = Ib_inv * (M + wbx*Ib*wb);

%% Accept delay using previous wb to calculate wi

euler_att   = EulerAngle_fromQuaternionData_NED(qb);
wi          = getInertialFromBody(wb,euler_att);

%% calc the differential quaternion

qb = qb/norm(qb);

dqb   = 0.5 * [ -wb'*qb(2:4); qb(1)*wb - wbx*qb(2:4)];


%% assemble the dx vector

dx = [dr;dv;dqb;dwb;de];


end