function dw = sim_body_vel(t,w,I_inv,I)

if(t<1)
    Mx = 100;
    My = 10;
    Mz = 1;
else
    Mx = 0;
    My = 0;
    Mz = 0;
end

M = [Mx;My;Mz];

wx = [ 0 -w(3) w(2);...
       w(3) 0 -w(1);...
       -w(2) w(1) 0];

dw = I_inv * (M + wx*I*w);

end