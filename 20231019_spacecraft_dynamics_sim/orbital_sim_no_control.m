
function [dx] = orbital_sim_no_control(t,x,ge,Ib,Ib_inv)

%% define variables
r = x(1:3);
qi = x(7:10);

%% first translation -- simple

dr = x(4:6);
dv = -ge * r * (r'*r)^-(3/2);

%% Then body rotation
% cannot find a nice solution for inertial rotation because the inertial
% tensor is not constant and would be a pain in the ass to calculate for
% each rotation

wb = x(11:13);

if(t<1)
    Mx = 0.1;
    My = 0.1;
    Mz = 0.1;
else
    Mx = 0;
    My = 0;
    Mz = 0;
end

M = [Mx;My;Mz];

wbx = [ 0 -wb(3) wb(2);...
       wb(3) 0 -wb(1);...
       -wb(2) wb(1) 0];

dwb = Ib_inv * (M + wbx*Ib*wb);

%% Accept delay using previous wb to calculate wi

euler_att   = EulerAngle_fromQuaternionData_NED(qi);
wi          = getInertialFromBody(wb,euler_att);

%% calc the differential quaternion

qi = qi/norm(qi);

dqi   = 0.5 * [ -wb'*qi(2:4); qi(1)*wb - wbx*qi(2:4)];

%% assemble the dx vector

dx = [dr;dv;dqi;dwb];


end