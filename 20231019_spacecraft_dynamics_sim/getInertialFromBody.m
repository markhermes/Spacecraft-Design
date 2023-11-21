function wi = getInertialFromBody(wb,euler)

    d2r = pi/180;

    r = d2r * euler(3);
    p = d2r * euler(2);
    y = d2r * euler(1);
    
    R = [1      0           -sin(p);...
        0       cos(r)      sin(r)*cos(p);...
        0       -sin(r)     cos(r)*cos(p)];

    if(size(wb,1) == 1)
        wi = flip(R*flip(wb'))';
    else 
        wi = flip(R*flip(wb))';
    end
end