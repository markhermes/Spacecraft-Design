function [err,tau] = get_fine_control_torques(Qfsn, Qb, Wb, Ib, DT, E_INT)

control_torques_max = 0.05;


for i = 1:length(Qfsn)

    wb      = Wb(i,:)';
    q_fsn   = Qfsn(i,:);
    qb      = Qb(i,:);

    wbx = [ 0 -wb(3) wb(2);...
        wb(3) 0 -wb(1);...
        -wb(2) wb(1) 0];

    % want the error in the body frame
    q_err = Quaternion_MULTIPLY(q_fsn,Quaternion_inv(qb));


    q_err = q_err/norm(q_err);


    % integrate the angular value for the torque correction... Direction should
    % always be applied in the error direction (assuming small angular
    % velocities)
    angular_error_mag = acos(abs(q_err(1)))/2 ;

    angular_error_dir = (q_err(2:4)'/norm(q_err(2:4)) - wb/norm(wb));

    e_val = angular_error_mag;
    e_int = E_INT(i);


    % maybe first just apply torques in the w_des direction?
    Kp = 0.0001;
    Ki = 0.000001;
    Kd = 0.001;
    control_torques = -(Kp*e_val + Ki*e_int)*angular_error_dir;

    % set max value of control_torques
    if(any(abs(control_torques) > control_torques_max))
        control_torques = control_torques_max* ...
            control_torques/max(abs(control_torques));
    end


    err(i,:) = 2*acos(abs(q_err(1)))*180/pi;
    tau(i,:) = control_torques';

end

end