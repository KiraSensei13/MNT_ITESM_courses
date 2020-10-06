function[x, y] = fsquarex(p3)
% p3 are the parameters of the square
% x0, y0 is the location of the center and a the apotema
% send parameters xo, yo and a
    n = 101;
    
    for i = 1 : n
        theta(i) = 2 * pi * (i - 1) / (n - 1);
        x0       = p3(1);
        y0       = p3(2);
        a        = p3(3);
        
        if theta(i) <= pi / 4 
            x(i) = x0 + a;
            y(i) = y0 + a * tan(theta(i));
        end
        
        if (theta(i) > pi / 4) && (theta(i) <= 3 * pi / 4)
            y(i) = y0 + a;
            x(i) = x0 + a * cot(theta(i));
        end
        
        if (theta(i) > 3 * pi / 4) && (theta(i) <= 5 * pi / 4)
            x(i) = x0 - a;
            y(i) = y0 + a * tan(pi - theta(i));
        end
        
        if (theta(i) > 5 * pi / 4) && (theta(i) <= 7 * pi / 4)
            y(i) = y0 - a;
            x(i) = x0 - a * cot(theta(i) - pi);
        end
        
        if theta(i) > 7 * pi / 4
            x(i) = x0 + a;
            y(i) = y0 + a * tan(theta(i));
        end
    end
    
    plot(x, y, 'ro-')
end