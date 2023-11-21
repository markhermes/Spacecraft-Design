%-------------------------------------------------------------------------
%
% 10-20-23 : Found that there is no reason to have the rotation in the
%               state matrix. Since the correct matrix is a body rate
%               transformation matrix, the error in the model will build
%               up. Also it is easy to confuse the rotation matrix for the
%               rate tranformation matrix, which appears to be the case
%               with current implementation in PHX
%
%       TODO: (1) use the realistic scenario with the accelerometers + heading
%       as the signal source. (2) try to fuse the ins rates with the gyro rates 
%
%-------------------------------------------------------------------------

%% pull in the required data for the filter to simulate
q_star_trk      = q_body;
ant_gyro_true   = y(:,11:13); %rad/s
time            = t;

DT  = t(10) - t(9);
d2r = pi/180;

%% noise the gyro measurements and show the uncorrected result

gyro_noise_mag = 0.001;

% add markov noise to the ant_gyro
x_noise = generateMarkovNoise(time);
y_noise = generateMarkovNoise(time);
z_noise = generateMarkovNoise(time);
ant_gyro = ant_gyro_true + [x_noise',y_noise',z_noise'] + gyro_noise_mag * randn(3,1);

% get Ant attitude from gyro only
q_ant_gyro_only_int = q_body(1,:);
for i = 2:numel(time)
    q_inc   = getQuatInc_fromGyroVectors(ant_gyro(i,:),DT);
    q_ant_gyro_only_int(i,:) =  Quaternion_MULTIPLY(q_inc,...
                                        q_ant_gyro_only_int(i-1,:));
    euler_gyro_only(i,:) = EulerAngle_fromQuaternionData_NED(...
                                        q_ant_gyro_only_int(i,:));
end

%% get magnetic field vectors (magnetometer)
% in reality would have to process the direction vector to a heading, but
% instead of reverse engineering that, just assume we can do it correctly.


heading_noise   = 0.1;

for i = 1:numel(time)

    euler_body_angle = EulerAngle_fromQuaternionData_NED(q_body(i,:));
    magnetic_yaw_inertial_to_bdy(i) = euler_body_angle(1) + heading_noise*randn;

end

%% setup the Kalman Filter matrices

gyro_correlation = -50;
gyro_bias_mat = gyro_correlation*eye(3);

% state transition model
F11 = 1*eye(3);
F21 = zeros(3);
F22 = 1*eye(3) + gyro_bias_mat;

% using e^AT = I + AT + (AT)^2/2
% F11 = 1*eye(3);
% F21 = zeros(3);
F12 = eye(3)*DT + gyro_bias_mat*DT^2/2;
F22 = eye(3) + gyro_bias_mat*DT + gyro_bias_mat^2*DT^2/2 + gyro_bias_mat^2*DT^3/6;
F   = [F11,F12;F21,F22];

% process noise model - gyros
Q = diag([ 0.001*ones(1,3) 0.05*ones(1,3)]);

% measurement model -- start with just INS
H_star_trk       = [eye(3) zeros(3)];
H_magnetometer   = [0 0 1 0 0 0]; %since the angles are


% measurement noise model
R_star_trk   = 0.001 * eye(3);
R_magnetometer   = 1;


% initialize the covariance matrix
P           = 1e-6 * eye(6);
state       = zeros(6,1);
gyro_bias_est = zeros(1,3);


%% run the kalman filter

magnetometer_update_period  = DT/DT;
star_trk_meas_update_period = floor(1/DT);

q_att = q_star_trk(1,:);
R_att = RotationMatrix_from_quaternion(q_att);

for i = 1:numel(time)

    euler_true_att(i,:) = EulerAngle_fromQuaternionData_NED(...
                                        q_body(i,:));
    
    % adjust the gyro bias based on the KF estimate
    ant_gyro_removed_bias(i,:) = ant_gyro(i,:) + (gyro_bias_est);
    
    % predict step is the integration of the gyros for 
    q_inc       = getQuatInc_fromGyroVectors(ant_gyro_removed_bias(i,:),DT);
    q_att       = Quaternion_MULTIPLY(q_inc,q_att);
    euler_att   = EulerAngle_fromQuaternionData_NED(q_att);
    
    % update covariance matrix
    P           = F*P*F' + Q;
    
    % run the correction step for startrk
    if(mod(i,star_trk_meas_update_period) == 0)
       
        q_err       = Quaternion_MULTIPLY(q_star_trk(i,:),Quaternion_inv(q_att));
        euler_err   = EulerAngle_fromQuaternionData_NED(q_err)
        
        % need to reverse euler err so it is roll, pitch, yaw so that the
        % gyro bias estimate will be multiplied correctly, since yaw = z,
        % pitch = y, roll = x
        euler_err   = flip(euler_err);
        
        H           = H_star_trk;
        K           = P*H'*(H*P*H' + R_star_trk)^-1;
        corrections = K*euler_err';
        state       = state + corrections
        P           = (eye(6)-K*H)*P;
        
        % correct attitude estimate
        q_correct = Quaternion_fromEulerAngle_YPR_NED(flip(corrections(1:3)));  
        
        % debug
        euler_ins   = EulerAngle_fromQuaternionData_NED(q_body(i,:))
        euler_att   = EulerAngle_fromQuaternionData_NED(q_att)
        q_att       = Quaternion_MULTIPLY(q_correct,q_att);
        euler_att  = EulerAngle_fromQuaternionData_NED(q_att)
        
        % est gyro bias by rotating the error --- effectively just a 1-d
        % kalman filter... not sure why the kalman filter method is not
        % working... Oh well this is fine...
        gyro_bias_est       = (state(4:6))';
        
        %disp('R mat');
        %disp(R_att);
        %disp('K mat');
        %disp(K(1:end,:));
        
    end 
    state(1:3) = zeros(3,1);

    %run the magnetometer correction -- quaternions to avoid euler shit
    if(mod(i,magnetometer_update_period) == 0)

        % create yaw only quaternion and meas quatern
        q_est_yaw = Quaternion_fromEulerAngle_YPR_NED([euler_att(1) 0 0]); 
        q_meas_yaw = Quaternion_fromEulerAngle_YPR_NED([euler_true_att(i,1) 0 0]); 
        % get the body error 
        q_err_yaw = Quaternion_MULTIPLY(q_meas_yaw,Quaternion_inv(q_est_yaw));

        % get the euler erro
        est_yaw   = EulerAngle_fromQuaternionData_NED(q_est_yaw);
        meas_yaw   = EulerAngle_fromQuaternionData_NED(q_meas_yaw);
        euler_err_yaw   = EulerAngle_fromQuaternionData_NED(q_err_yaw);
        euler_err       = [euler_err_yaw(1)]

        H           = H_magnetometer;
        K           = P*H'*(H*P*H' + R_magnetometer)^-1;
        corrections = K*euler_err';
        state       = state + corrections
        P           = (eye(6)-K*H)*P;

        % correct attitude estimate
        yaw_correct = corrections(3);
        euler_att   = EulerAngle_fromQuaternionData_NED(q_att)
        euler_att(1) = euler_att(1) + yaw_correct;
        q_att       = Quaternion_fromEulerAngle_YPR_NED(euler_att);

        % debug
        euler_ins       = EulerAngle_fromQuaternionData_NED(q_body(i,:))
        corrected_att   = EulerAngle_fromQuaternionData_NED(q_att)

        % calculate the error
        euler_est_att(i,:) = EulerAngle_fromQuaternionData_NED(q_att);
        kf_attitude_error(i,:) = norm(euler_est_att(i,:) - euler_true_att(i,:));


        % est gyro bias by rotating the error --- effectively just a 1-d
        % kalman filter... not sure why the kalman filter method is not
        % working... Oh well this is fine...
        gyro_bias_est       = (state(4:6))';

        %disp('R mat');
        %disp(R_att);
        %disp('K mat');
        %disp(K(1:end,:));

    end    
    state(1:3) = zeros(3,1);

    
end


uncorrected_euler_err = vecnorm( (euler_true_att -  euler_gyro_only)' ) ;

plot_kf_performance;


%% FUNCTIONS
function R = RotationInertial_fromBodyRates(angles)

y = d2r * angles(1); p = d2r * angles(2); r = d2r * angles(3);

R = [1 sin(r)*tan(p) cos(r)*tan(p);...
    0 cos(r) -sin(r);...
    0 sin(r)*sec(p) cos(r)*sec(p)];

end

function R = RotationBody_fromInertialRates(angles)

y = d2r * angles(1); p = d2r * angles(2); r = d2r * angles(3);

R = [1 0 0;...
    0 cos(r) sin(r)*cos(p);...
    0 -sin(r) cos(r)*cos(p)];


end

function q = Quaternion_fromAxisAngle(axis_angle)

w = sqrt(1-axis_angle'*axis_angle);
q = [w axis_angle'];

end

% accepts rad/s gyro vectors
function q = getQuatInc_fromGyroVectors(v,DT)

    mag = norm(v);
    angle = mag*DT*0.5;
    
    q(1) = cos(angle);
    q(2:4) = v * sin(angle)/mag;

end

function v = getGyroVectors_fromQuatInc(q,DT)

    angle   = acos(q(1));
    mag     = angle/DT*2; 
    v       = q(2:4)*mag/sin(angle);
    
end


function noise = generateMarkovNoise(t)

% Parameters
n = numel(t); % Length of the Markov noise sequence
mu = 0.005; % Mean of the noise
sigma = 0.00001; % Standard deviation of the noise
rho = 1; % Autocorrelation coefficient (between -1 and 1)

% Generate Markov noise
noise = zeros(1, n); % Initialize the noise vector

% Generate the first value from a normal distribution
noise(1) = mu + sigma * randn(1);

% Generate the rest of the values based on the Markov process
for i = 2:n
    noise(i) = mu + rho * (noise(i-1) - mu) + sigma * randn(1);
end

% Plot the Markov noise
% plot(noise)
% % title('Markov Noise')
% xlabel('Time')
% ylabel('Value')

end
