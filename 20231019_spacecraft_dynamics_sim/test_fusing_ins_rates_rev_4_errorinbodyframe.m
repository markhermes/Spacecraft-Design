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

clear all;

d2r = pi/180;

motion_freq         = 0.01;
pitch_flutter_freq  = 10;
roll_flutter_freq   = 20;
pitch_offset        = 0.1;
DT                  = 0.001;
time                = 0:DT:100;

%% create the ins motion and the quaternion increment
yaw     = 120*sin(2*pi*motion_freq*time);
pitch   = 30*sin(2*pi*motion_freq*time);
roll    = 10*sin(pi*motion_freq*time);
q_ins_int_test  = [1 0 0 0];
for i = 1:numel(time)
    q_ins(i,:) = Quaternion_fromEulerAngle_YPR_NED([yaw(i) pitch(i) roll(i)]);
    
    %calc backward deriv increment
    if(i > 1)
        q_ins_inc(i,:) = Quaternion_MULTIPLY(q_ins(i,:),...
            Quaternion_inv(q_ins(i-1,:)));
                
        %make sure the result is the same as q_ant
        q_ins_int_test(i,:) = Quaternion_MULTIPLY(q_ins_inc(i,:),q_ins_int_test(i-1,:));
    end
end

% get INS gyro vectors
for i = 2:numel(time)
    mag         = acos(q_ins_inc(i,1))/(0.5*DT);
    half_angle  = mag*DT*0.5;
    ins_gyro(i,:)         = q_ins_inc(i,2:4)*mag/sin(half_angle);
end


%% create the bending offsets
% flutter is occuring in the local level frame... so first create the
% motion quat then add on the flutter...

bending_yaw     = 0*sin(2*pi*time);
bending_pitch   = pitch_offset + 0.5*sin(2*pi*pitch_flutter_freq*time);
bending_roll    = 0.1*sin(2*pi*roll_flutter_freq*time);
q_ant_int_test  = [1 0 0 0];
for i = 1:numel(time)
    q_wing      = Quaternion_fromEulerAngle_YPR_NED(...
        [bending_yaw(i) bending_pitch(i) bending_roll(i)]);
    
    q_ant(i,:)  = Quaternion_MULTIPLY(q_wing, q_ins(i,:));
    euler_ant(i,:) = EulerAngle_fromQuaternionData_NED(q_ant(i,:));
    
    if(i > 1)
       q_ant_inc(i,:) = Quaternion_MULTIPLY(q_ant(i,:),...
            Quaternion_inv(q_ant(i-1,:)));
        
        %make sure the result is the same as q_ant
        q_ant_int_test(i,:) = Quaternion_MULTIPLY(q_ant_inc(i,:),q_ant_int_test(i-1,:));
    end
end

%% get Ant gyro vectors - euler then convert to body

for i = 2:numel(time)
    mag         = acos(q_ant_inc(i,1))/(0.5*DT);
    half_angle  = mag*DT*0.5;
    ant_gyro_true(i,:)         = q_ant_inc(i,2:4)*mag/sin(half_angle);
end

% add markov noise to the ant_gyro
x_noise = generateMarkovNoise(time);
y_noise = generateMarkovNoise(time);
z_noise = generateMarkovNoise(time);
ant_gyro = ant_gyro_true + [x_noise',y_noise',z_noise'];

% get Ant attitude from gyro only
q_ant_gyro_only_int = [1 0 0 0];
for i = 2:numel(time)
    q_inc   = getQuatInc_fromGyroVectors(ant_gyro(i,:),DT);
    q_ant_gyro_only_int(i,:) =  Quaternion_MULTIPLY(q_inc,...
                                        q_ant_gyro_only_int(i-1,:));
    euler_gyro_only(i,:) = EulerAngle_fromQuaternionData_NED(...
                                        q_ant_gyro_only_int(i,:));
                                    
end

%% get Down Vectors (accelerometer)

% simulate a jet rapid acceleration / decceleration, idk...
acc_noise_freq  = 0.01;
acc_noise_mag   = 0.02;

for i = 1:numel(time)
    v_down(i,:) = Quaternion_rotateVectorWithQuat([0 0 0 1],q_ant(i,:));
    
    % add noise... maybe sinusoidal? -- based on acceleration ...
    v_down(i,1) = v_down(i,1) + acc_noise_mag*sin(acc_noise_freq*2*pi*time(i));
    v_down(i,2) = v_down(i,2) + acc_noise_mag*sin(acc_noise_freq*2*pi*time(i));
    v_down(i,3) = v_down(i,3) + acc_noise_mag*sin(acc_noise_freq*2*pi*time(i));

end

%% get yaw (ins)
heading = euler_ant(:,1);

%% setup the Kalman Filter matrices

gyro_correlation = -100;
gyro_bias_mat = gyro_correlation*eye(3);

% state transition model
F11 = 1*eye(3);
F21 = zeros(3);
F22 = 1*eye(3) + gyro_bias_mat;

% using e^AT = I + AT + (AT)^2/2
F11 = 1*eye(3);
F21 = zeros(3);
F12 = eye(3)*DT + gyro_bias_mat*DT^2/2;
F22 = eye(3) + gyro_bias_mat*DT + gyro_bias_mat^2*DT^2/2;
F   = [F11,F12;F21,F22];

% process noise model - gyros
Q = diag([ 0.001*ones(1,3) 0.01*ones(1,3)]);

% measurement model -- start with just INS
H   = [eye(3) zeros(3)];

% measurement noise model
R   = 0.001 * eye(3);

% initialize the covariance matrix
P           = 0.1 * eye(6);
state       = zeros(6,1);
gyro_bias_est = zeros(1,3);


%% run the kalman filter
state = [0;0;0;0;0;0];
q_att = q_ant(1,:);
R_att = RotationMatrix_from_quaternion(q_att);
meas_update_period = 40;
total_err = 0;
for i = 1:numel(time)

    % adjust the gyro bias based on the KF estimate
    ant_gyro_removed_bias(i,:) = ant_gyro(i,:) + (gyro_bias_est);
    
    % predict step is the integration of the gyros for 
    q_inc       = getQuatInc_fromGyroVectors(ant_gyro_removed_bias(i,:),DT);
    q_att       = Quaternion_MULTIPLY(q_inc,q_att);
    euler_att   = EulerAngle_fromQuaternionData_NED(q_att);
    
    % update covariance matrix
    P           = F*P*F' + Q;
    
    % run the correction step
    if(mod(i,meas_update_period) == 0)
       
        q_err       = Quaternion_MULTIPLY(q_ant(i,:),Quaternion_inv(q_att));
        euler_ins   = EulerAngle_fromQuaternionData_NED(q_ant(i,:))
        euler_att   = EulerAngle_fromQuaternionData_NED(q_att)
        euler_err   = EulerAngle_fromQuaternionData_NED(q_err)
        
        % need to reverse euler err so it is roll, pitch, yaw so that the
        % gyro bias estimate will be multiplied correctly, since yaw = z,
        % pitch = y, roll = x
        euler_err   = flip(euler_err);
        
        K           = P*H'*(H*P*H' + R)^-1;
        corrections = K*euler_err';
        state       = state + corrections
        P           = (eye(6)-K*H)*P;
        
        % correct attitude estimate
        q_correct = Quaternion_fromEulerAngle_YPR_NED(flip(corrections(1:3)));  
        
        % debug
        q_att       = Quaternion_MULTIPLY(q_correct,q_att);
        euler_test  = EulerAngle_fromQuaternionData_NED(q_att)
        
        % est gyro bias by rotating the error --- effectively just a 1-d
        % kalman filter... not sure why the kalman filter method is not
        % working... Oh well this is fine...
        gyro_bias_est       = (state(4:6))';

        total_err = total_err + norm(state(1:3));
        
        disp('R mat');
        disp(R_att);
        disp('K mat');
        disp(K(1:end,:));
        
    end 
end
    for j = 1:3  subplot(1,3,j);  plot(time,ant_gyro_removed_bias(:,j),'.k'); hold on;
    plot(time,ant_gyro_true(:,j),'.r'); end
    xlim([0 100]);
    title(sprintf('total error %0.1f',total_err))


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
mu = 0.2; % Mean of the noise
sigma = 0.001; % Standard deviation of the noise
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
