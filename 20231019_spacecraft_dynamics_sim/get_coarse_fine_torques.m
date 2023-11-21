function [err,tau] = get_fine_control_torques(Qfsn, Qb, Wb, Ib, DT)

control_torques_max = 0.1; 


for i = 1:length(Qfsn)

    wb      = Wb(i,:)';
    q_fsn   = Qfsn(i,:);
    qb      = Qb(i,:);

    wbx = [ 0 -wb(3) wb(2);...
        wb(3) 0 -wb(1);...
        -wb(2) wb(1) 0];

    % want the error in the body frame
    q_err = Quaternion_MULTIPLY(q_fsn,Quaternion_inv(qb));

    % the axis of the quaternion defines the optimal w vector
    if(q_err(1) < 0) q_err = -q_err; end
    w_des_dir = q_err(2:4)';
    q_err(1);

    % scale the desired vector by how large the error is
    max_w_des = 0.2;
    w_des = w_des_dir * max_w_des * (1.5-abs(q_err(1)));

    % maybe first just apply torques in the w_des direction?
    control_torques = (Ib*(w_des - wb)/DT - wbx * Ib * wb);

    % set max value of control_torques
    if(any(abs(control_torques) > control_torques_max))
        control_torques = control_torques_max* ...
            control_torques/max(abs(control_torques));
    end


    err(i,:) = 2*acos(q_err(1))*180/pi;
    tau(i,:) = control_torques';

end

end