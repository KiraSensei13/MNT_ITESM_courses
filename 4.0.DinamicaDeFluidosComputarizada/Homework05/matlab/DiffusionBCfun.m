function [pl, ql, pr, qr] = DiffusionBCfun(xl, u1, xr, ur, t, P)
    % Boundary conditions for x = 0 and x = L;
    % pL(xl, t, ul) + qL(xl, t) * f(x, t, u, du / dx) = 0 % Left boundary
    % pR(xr, t, ur) + qR(xr, t) * f(x, t, u, du / dx) = 0 % Right boundary
    
    % Extract parameters
    c0 = P(2);
    
    % BCs: No flux boundary at the right boundary and constant
    % concentration on the left boundary
    pl = u1 - c0;
    ql = 0;
    pr = 0;
    qr = 1;
end