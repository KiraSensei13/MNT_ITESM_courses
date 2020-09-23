function[ydot] = drainingTank(x, y, par)
    g      = par(1);
    lambda = par(2);
    L      = par(3);
    
    ydot(1) = y(2);
    ydot(2) = (-g*y(1) + ((y(2))^2*(lambda^4 - 1))/2) / (L*lambda^2  + y(1));
    ydot    = ydot';
end