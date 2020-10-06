function[xi,time,Y]=ReactorPFDD() 
p(1) = 0.0005;   %Diffusion coefficient D
 p(2) = 1;      % c0
 p(3)=1.5;      % k , first order kinetic coefficient
 p(4)=1;        % vo , velocity of fluid injection
 M=60;
 % This can be nn=2,3 or 4
 nn=3;
 p(5)=M;
 p(6)=nn;        % Truncation error power
Tspan=[0 1];
tspan=linspace(0,1,20);
xi=linspace(0,1,M);
Y0=zeros(M,1);
Y0(1)=1.0;
OPTIONS=[];
[time,Y]=ode45(@reactub1,tspan,Y0,[],p);
subplot(2,1,1);
plot(time,Y);
xlabel('\tau');ylabel('Concentration [mol/dm^3]')
subplot(2,1,2);
plot(xi,Y');
xlabel('distance x/L');ylabel('Concentration [mol/dm^3]')
end
% tubular reactor by jose lopez salinas
function yprime=reactub1(t,y,p)
c=y;
 D = p(1);
 k=p(3);
 vo=p(4);
 N=p(5);
 nn=p(6);
 m=1;
yprime(1)=0;
dx=1/(N-1);
if nn==2
    for i=2:N-1
    yprime(i)=D*(c(i+1)-2*c(i)+c(i-1))/(dx^2)-vo*(c(i+1)-c(i-1))/(2*dx)-k*c(i)^m;
    end
    yprime(N)=(4*yprime(N-1)-yprime(N-2))/3;
elseif nn==3
    % ====================== o (h^4)=====================
     for i=3:N-2
       yprime(i)=D*(-c(i+2)/12+4*c(i+1)/3-5*c(i)/2+4*c(i-1)/3-c(i-2)/12)/(dx^2);
       yprime(i)=yprime(i)-vo*(c(i-2)/12-2*c(i-1)/3+2*c(i+1)/3-c(i+2)/12)/dx;
       yprime(i)=yprime(i)-k*c(i)^m;
     end
      yprime(2)=D*(11*c(1)-20*c(2)+6*c(3)+4*c(4)-c(5))/(12*dx^2);
      yprime(2)=yprime(2)-vo*(-2*c(1)-3*c(2)+6*c(3)-c(4))/(16*dx);
      yprime(2)=yprime(2)-k*c(2)^m;

      yprime(N-1)=D*(11*c(N)-20*c(N-1)+6*c(N-2)/6+4*c(N-3)-c(N-4))/(12*dx^2);
      yprime(N-1)=yprime(N-1)-vo*(c(N-3)-6*c(N-2)+3*c(N-1)+2*c(N))/(6*dx);
      yprime(N-1)=yprime(N-1)-k*c(N-1)^m;
      % Boundary condition
      yprime(N)=(18*yprime(N-1)-9*yprime(N-2)+2*yprime(N-3))/11;
      %================================================
else
% ====================== o (h^4)=====================
     for i=3:N-2
       yprime(i)=D*(-c(i+2)/12+4*c(i+1)/3-5*c(i)/2+4*c(i-1)/3-c(i-2)/12)/(dx^2);
       yprime(i)=yprime(i)-vo*(c(i-2)/12-2*c(i-1)/3+2*c(i+1)/3-c(i+2)/12)/dx;
       yprime(i)=yprime(i)-k*c(i)^m;
     end
      yprime(2)=D*(15*c(2)/4-77*c(3)/6+107*c(4)/6-13*c(5)+61*c(6)/12-5*c(7)/6)/(dx^2);
      yprime(2)=D*(10*c(1)-15*c(2)-4*c(3)+14*c(4)-6*c(5)+c(6))/(12*dx^2);
      %yprime(2)=yprime(2)-vo*(-25*c(2)/12+4*c(3)-3*c(4)+4*c(5)/3-c(6)/4)/dx;
      yprime(2)=yprime(2)-vo*(-3*c(1)-10*c(2)+18*c(3)-6*c(4)+c(5))/(12*dx);
      yprime(2)=yprime(2)-k*c(2)^m;

      yprime(N-1)=D*(15*c(N-1)/4-77*c(N-2)/6+107*c(N-3)/6-13*c(N-4)+61*c(N-5)/12-5*c(N-6)/6)/(dx^2);
      yprime(N-1)=D*(c(N-5)-6*c(N-4)+14*c(N-3)-4*c(N-2)-15*c(N-1)+10*c(N))/(12*dx^2);
      %yprime(N-1)=yprime(N-1)-vo*(25*c(N-1)/12-4*c(N-2)+3*c(N-3)-4*c(N-4)/3+c(N-5)/4)/dx;
      yprime(N-1)=yprime(N-1)-vo*(-c(N-4)+6*c(N-3)-18*c(N-2)+10*c(N-1)+3*c(N))/(12*dx);
      yprime(N-1)=yprime(N-1)-k*c(N-1)^m;
      % Boundary condition
      yprime(N)=4*yprime(N-1)-3*yprime(N-2)+4*yprime(N-3)/3-yprime(N-4)/4;
      yprime(N)=12*yprime(N)/25;
      %yprime(N)=15*yprime(N-1)-6*yprime(N-2)+yprime(N-3);
      %yprime(N)=yprime(N)/10;
      %================================================
 end
 yprime=yprime';
end
