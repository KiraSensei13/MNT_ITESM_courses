function [c, f, s] = DiffusionPDEfun(x, t, u, dudx, P)
    % Function defining the PDE
    % C du/dt =(1/x^n)d(x^n f)/dx+s
    % f=f(x,t,u,du/dx)
    % s=s(x,t,u,du/dx)
    % c=c(x,t,u,du/dx)
    
    % Extract parameters
    D  = P(1);
    k  = P(3);
    vo = P(4);
    
    % PDE
    c = 1;
    f = D .* dudx;
    s = - k * u - vo * dudx;
end