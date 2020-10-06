% Function to calcualte velocity of each parcel
function[zprime] = Vfield(t, z, p)
    [n1,m1] = size(z);
    n       = max(n1, m1);
    
    for i = 1 : n / 2
        x(i) = z(2 * i - 1);
        y(i) = z(2 * i);
    end
    
    for i = 1 : n / 2
        Vx(i) = 2 * x(i) * y(i);
        Vy(i) = - y(i) ^ 2;
    end
    
    for i = 1 : n / 2
        zprime(2 * i - 1) = Vx(i);
        zprime(2 * i)     = Vy(i);
    end
    
    zprime = zprime';
end