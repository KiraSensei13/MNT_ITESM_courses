function [xx, yy] = lineFunction(Ax, Ay, Bx, By, n)
    xx        = linspace(Ax, Bx, n);
    slope     = (By - Ay)/(Bx - Ax);
    intercept = Ay - slope * Ax;
    yy        = slope * xx + intercept;
end

