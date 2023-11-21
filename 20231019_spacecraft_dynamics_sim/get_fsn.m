function q_fsn = get_fsn(R,V)

for i = 1:length(R)

    r = R(i,:);
    v = V(i,:);

    % forward is velocity vector
    F = v/norm(v);

    % nadir is -r
    N = -r/norm(r);

    % starboard is -FxN
    S = -cross(F,N);

    R_fsn = [F',S',N'];
    q_fsn(i,:) = Quaternion_fromRotationMatrix(R_fsn);

end

end