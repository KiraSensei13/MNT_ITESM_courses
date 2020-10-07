% function to draw a circle
function[x, y] = fcircle(p1)
    n     = 40;
    m     = n + 1;
    beta1 = 2 * pi / 1;
    x0    = p1(1);
    y0    = p1(2);
    rho   = p1(3);
    
    % use parameters x0, y0 and rho (center of the circle and radius)
    for i = 1 : m
        theta(i) = (i-1) * beta1 / n;
        x(i)     = rho * cos(theta(i));
        y(i)     = rho * sin(theta(i));
    end
    
    for i = 1 : m
        x(i) = x0 + x(i);
        y(i) = y0 + y(i);
    end
end