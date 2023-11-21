function [value, isterminal, direction] = turn_on_fine_control(t,y)

    value = norm(y(11:13)) > 0.01;

end