% Program to se the velocity field and to track a tracer figure
% By JLLopez CFD 09/29/2020
% The tracer figure is a circle
rho=0.125;              % Raduius of the circle
 x0=0.95;y0=0.75;     % Center of the circle
 p1=[x0 y0 rho];      % parameter to draw the circle
%[x,y]=fsquarex(p1); % function to generate the shape
 [x,y]=fcircle(p1);      % function to generate the shape
    [n1,m1]=size(x);
    m=max(n1,m1);
    to=0;tf=0.100;
    % vector position
    for i=1:m
        z(2*i-1)=x(i);
        z(2*i)=y(i);
    end
    % parameters you may need in the vector field function
    p=1;
    zo=z;  % initial condition
    tspan=linspace(0,0.25,20); % time span to track the fluid parcels
    tf=[0 1];
    % Solution of the ODEs dr/dt , here you solve the velocity field eqn (vector field) 
    [time,YS]=ode45(@Vfield,tspan,zo,[],p);
    % You will plot that last time position
    YL=YS(end,:);
    for i=1:m
        xf(i)=YL(2*i-1);
        yf(i)=YL(2*i);
    end
    [mt,nt]=size(time);
    mt=max(nt,mt);
    for j=1:mt
        for i=1:m
            XP(j,i)=YS(j,2*i-1);
            YP(j,i)=YS(j,2*i);
        end
    end
    % p3=p1;
    % [xs,ys]=fsquarex(p3);
    %plot(xs,ys,'rs-');
    %hold;
    % plots arrows with directional components U and V at the Cartesian coordinates specified by X and Y
    hold on;
    [x2,y2,Ux,Uy]=MPlotxx1();
    % plots initial position and final to compare
    plot(x,y,'b.-',xf,yf,'r.-');
    xlabel('x','fontsize',18);ylabel('y','fontsize',18);
    xlim([0.0 2.0]);ylim([0.0 1.5]);
    colorbar;
% function to draw a circle
for i=2:mt
    title(sprintf('time=%1.4f',time(i)),'fontsize',16);
    plot(XP(i,:),YP(i,:),'r.-');
    shg;
    pause(1);
end
function[x,y]=fcircle(p1)
    n=40;
    m=n+1;
    beta1=2*pi/1;
    x0=p1(1);y0=p1(2);rho=p1(3);
    % use parameters x0, y0 and rho (center of the circle and radius)
    for i=1:m
        theta(i)=(i-1)*beta1/n;
        x(i)=rho*cos(theta(i));
        y(i)=rho*sin(theta(i));
    end
    for i=1:m
        x(i)=x0+x(i);
        y(i)=y0+y(i);
    end
end
function[x,y]=fsquarex(p3)
% p3 are the parameters of the square
% x0, y0 is the location of the center and a the apotema
% send parameters xo, yo and a
    n=101;
    for i=1:n
        theta(i)=2*pi*(i-1)/(n-1);
    x0=p3(1);y0=p3(2);a=p3(3);
    if theta(i)<=pi/4 
        x(i)=x0+a;
        y(i)=y0+a*tan(theta(i));
    end
    if (theta(i)>pi/4)&&(theta(i)<=3*pi/4)
        y(i)=y0+a;
        x(i)=x0+a*cot(theta(i));
    end   
    if (theta(i)>3*pi/4)&&(theta(i)<=5*pi/4) 
        x(i)=x0-a;
        y(i)=y0+a*tan(pi-theta(i));
    end
    if (theta(i)>5*pi/4)&&(theta(i)<=7*pi/4) 
        y(i)=y0-a;
        x(i)=x0-a*cot(theta(i)-pi);
    end 
        if theta(i)>7*pi/4 
        x(i)=x0+a;
        y(i)=y0+a*tan(theta(i));
    end
end
%plot(x,y,'ro-')
end
% Function of the vector field
function[X,Y,Ux,Uy]=MPlotxx1()
   [X,Y] = meshgrid(0.01:0.05:2.0,0.01:0.05:2.0);
     % Ux=2*X.*Y;Uy=-Y.^2;
     Ux = +X; %X.^2-Y.^2;
     Uy = +Y; %-1-2*X.*Y;
     quiver(X,Y,Ux,Uy,2);
     %phi=(Y.^3)/3-(X.^2).*Y;
     phi=-(X.^3)/3+X.*Y.^2+Y;
     contour(X,Y,phi);
     % [m,n]=size(X);
     % for i=1:m
      %  for j=1:n
      %       phi(i,j)=Y(i,j)^3/3-X(i,j)^2*Y(i,j);
      %   end
      % end
      % contour(X,Y,phi);
end
% Function to calcualte velocity of each parcel
function[zprime]=Vfield(t,z,p)
    [n1,m1]=size(z);
    n=max(n1,m1);
    for i=1:n/2
        x(i)=z(2*i-1);
        y(i)=z(2*i);
    end
    for i=1:n/2
         %Vx(i)=2*x(i)*y(i);
         %Vy(i)=-y(i)^2;
         Vx(i)= +x(i); %x(i)^2-y(i)^2;
         Vy(i)= +y(i); %-1-2*x(i)*y(i);
    end
    for i=1:n/2
        zprime(2*i-1)=Vx(i);
        zprime(2*i)=Vy(i);
    end
    zprime=zprime';
end