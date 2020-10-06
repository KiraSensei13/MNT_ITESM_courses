function yprime=TankDENV(t,y,p)
% This is the transicent of inviscid flow within a tank
% This just to make sure we are dealing with column vector
yprime=[0 0]';
%p=[DT,d,L,g,b];
DT=p(1);     % [m], diameter of the tank
d=p(2);   %[m], diameter of the pipe
L=p(3);         %[m], length of the pipe
g=p(4);       %[m/s^2], gravity 
b=p(5);     % Ratio of the diameters
% In in this case height is the dependent variable y1=h
yprime(1)=y(2);
yprime(2)=(-g*y(1)+((b^4-1)*y(2)^2)/2)/(b^2*L+y(1));
if y(1)<0
    yprime(1)=0;
    yprime(2)=0;
end

