function[ydot] = growingBubbles(x, y, par)
    k   = par(1);
    a1  = par(2);
    b1  = par(3);
    tau = x;
    B   = (1+cos(pi*tau/5))/2;
    if tau > 10
        B = 1;
    end
    
    ydot(1) = y(2);
    ydot(2) = (-B+1/(y(1))^(3*k)-(3/2)*y(2)^2-a1*y(2)/y(1)-b1/y(1))/y(1);
    ydot    = ydot';
end