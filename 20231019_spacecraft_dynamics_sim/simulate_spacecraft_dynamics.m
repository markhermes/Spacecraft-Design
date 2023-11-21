% 10/19 - todo : derive the inertia tensor using rudimentary equations.
% 10/20 - todo : see what happens when inertia matrix is changed.

%% reset data

close all; clear all;

%% run simulation for angular velocity

% assume translation is zero
trans = [0,0,0];

% create inertia matrix for cylinder
Izz = 0.1;
Ixx = 10;
Iyy = 10;
I = diag([Ixx, Iyy, Izz]);
I_inv = inv(I);

% small torque about the z-axis
% equation is wbdot = Iinv ( torque - wb x I x wb)
wb0 = [0;0;0];

opts = odeset('MaxStep',1e-3);
[t,wb] = ode45(@(t,w) sim_body_vel(t,w,I_inv,I),[0 5],[0,0,0],opts);

q_state = Quaternion_fromEulerAngle_YPR_NED([0 45 0]);
for i = 2:length(wb)

    % first get body rates into inertial rates...
    euler_att   = EulerAngle_fromQuaternionData_NED(q_state);
    wi          = getInertialFromBody(wb(i,:),euler_att);
    wi_vect(i,:)= wi';

    % next integrate rates
    DT          = t(i)-t(i-1);
    dq          = getQuatInc_fromGyroVectors(wb(i,:),DT);
    q_state     = Quaternion_MULTIPLY(dq,q_state);
    q_body(i,:) = q_state;
    euler_angles(i,:) = euler_att;
end


%% run movie

for i = 2:100:length(wb)

    % plot spacecraft
    plot_space_craft(q_body(i,:),trans);
    pause(0.01);
end

%% functions 

function wi = getInertialFromBody(wb,euler)

    d2r = pi/180;

    r = d2r * euler(3);
    p = d2r * euler(2);
    y = d2r * euler(1);
    
    R = [1      0           -sin(p);...
        0       cos(r)      sin(r)*cos(p);...
        0       -sin(r)     cos(r)*cos(p)];

    wi = flip(R*flip(wb'))';

end



function q = getQuatInc_fromGyroVectors(v,DT)

    mag = norm(v);
    angle = mag*DT*0.5;
    
    q(1) = cos(angle);
    q(2:4) = v * sin(angle)/mag;

end