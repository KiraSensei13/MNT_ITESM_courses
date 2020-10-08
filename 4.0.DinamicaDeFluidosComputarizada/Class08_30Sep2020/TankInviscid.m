Tspan=[0 1];
tspan=linspace(0,5.5,150);
Y0(1)=1;
Y0(2)=0;
DT=0.9095;     % [m], diameter of the tank
%DT=0.5;     % [m], diameter of the tank
d=2.5e-2;   %[m], diameter of the pipe
L=1.0;         %[m], length of the pipe
g=9.81;       %[m/s^2], gravity 
b=DT/d;     % Ratio of the diameters
p=[DT,d,L,g,b];
[time,Y]=ode45(@TankDENV,tspan,Y0,[],p);
plot(time,Y(:,1),time,-b^2*Y(:,2))
xlabel('time(s)')
ylabel('height')
title('Tank height and velocity')
