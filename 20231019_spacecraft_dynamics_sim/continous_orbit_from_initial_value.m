% Medium satellite	1,201 to 2,500 kg
% Small satellite	601 to 1,200 kg
% 
try close(v);
catch
end

clear all;

%% Create inertia matrix for cylinder -- dimensions of ISS 
% solid cylinder - Izz = 1/2 mr^2, Ixx = Iyy = 1/12(3r^2 + h^2)
% assume a medium satellite has dims of 1m by 0.2m
sat_m = 2e2;
sat_r = 0.2;
sat_l = 1;
Izz = 0.5 * sat_m * sat_r^2;
Ixx = 1/12 * sat_m * (3*sat_r^2 + sat_l^2);
Iyy = Ixx;
Ib = diag([Ixx, Iyy, Izz]);
Ib_inv = inv(Ib);

%% Get orbital path

% F = G m1 m2/r^2 , G m1 = ge, F = ge m2 / r^2
ge = 398600.8; % Earth gravitational constant [km3/s2]

%
orbital_body_mass   = 420e3;        %kg
earth_mass          = 5.9722e24;    %kg

% two inertial coord frames from TLE - ECI, TEME - True Equator Mean
% Equinox -- z axis is celestial ephemeris pole
% r is km, v = km/s
r_orbital_body   = [5.097299027407659e+03;3.155106125548773e+03;-3.209936007009090e+03]; %r_teme
v_orbital_body   = [-0.655709206595218;5.930726505531845;4.795677215184148]; %v_teme

% simulate the angular and translational dynamics
q0      = Quaternion_fromEulerAngle_YPR_NED([0, 45 , 0])';
wb0     = [0 0 0];   % initial z- spin
tspan   = [0 1000]; % orbital period is 90 mins
ic      = [r_orbital_body ; v_orbital_body ; q0'; wb0'];
opts    = odeset('MaxStep',0.1);
[t,y]   = ode45(@(t,x) orbital_sim_no_control(t,x,ge,Ib,Ib_inv),tspan,ic,opts);

r_orbital_body  = y(:,1:3);
v_orbital_body  = y(:,4:6);
q_body          = y(:,7:10);

% show convergence using maxstep of 0.1 secs
% display(y(end,:)')

% normalize the translation values for the sim
trans   = r_orbital_body/(1.5*max(max(r_orbital_body)));
vel     = v_orbital_body/(1.5*max(max(r_orbital_body)));

kf_body_frame;


%% Run Kalman Filter



%% run movie
v = VideoWriter('orbit_with_control.mp4','MPEG-4');
v.set('FrameRate',30,'Quality',100)
open(v);

for i = 2:(10*60):length(t)



    % plot spacecraft
    plot_space_craft(q_body(i,:),trans(i,:));

    % plot orbit
    plot3(trans(:,1),trans(:,2),trans(:,3),'-r'); hold on;

    % plot primary body
    plot_earth;

    % plot ascending and descending nodes
    plot_nodes(trans);

    % plot ForwardStarboardNadir
    plot_fsn(trans(i,:),vel(i,:));

    pause(0.01);

    writeVideo(v,getframe(gcf));


end
close(v);
