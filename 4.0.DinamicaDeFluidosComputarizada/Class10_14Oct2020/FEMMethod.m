ao=[1 1];
PHI=3.0;
s=PHI;
x1=1/3;
x2=2/3;
% x1=0.25;
% x2=0.75;
p=[PHI x1 x2];
[a,fval]=fsolve(@ftsolve1,ao,[],p);
NP=20;
xd=linspace(0.01,0.99,NP);
yan=sinh(xd*PHI)./(xd*sinh(PHI));
% plot(xd,yan);
for i=1:NP
    x1=xd(i);
    [TF,TFD,TFDD]=TestFuns(x1);
    y1=TF(1)+a(1)*TF(2)+a(2)*TF(3);
    yapp(i)=y1;
end
title(sprintf('PHI=%1.4f',s,'fontsize',16));
plot(xd,yan,xd,yapp,'ro');
xlabel('r/R');
ylabel('C/Cs');
title(sprintf('PHI=%1.4f',PHI));
%title(sprintf('PHI=%1.4f',s,'fontsize',16));
%======================================
function [F]=ftsolve1(a,p)
x1=p(2);x2=p(3);
phi=p(1);
x(1)=x1;
x(2)=x2;
%===================
% f =phi(1)+a(1)*phi(2)+a(3)*phi(3)
% DETS: x^2*d(dy/dx)/dx+2*x*dy/dx-PHI^2*x^2*y =0
for i=1:2 % two internal points
    xd=x(i);
    [g,gd,gdd]=TestFuns(xd);
    for j=1:3 % three test functions
    TF(i,j)=g(j);
    TFD(i,j)=gd(j);
    TFDD(i,j)=gdd(j);
    end
    fun(i)=TF(i,1);
    fund(i)=TFD(i,1);
    fund2(i)=TFDD(i,1);
    for j=2:3
        fun(i)=fun(i)+a(j-1)*TF(i,j);
        fund(i)=fund(i)+a(j-1)*TFD(i,j);
        fund2(i)=fund2(i)+a(j-1)*TFDD(i,j);
    end
    % x^2*d(dy/dx)/dx+2*x*dy/dx-PHI^2*x^2*y =0
    RES(i)=x(i)^2*fund2(i)+2*x(i)*fund(i)-phi^2*x(i)^2*fun(i);
    F(i)=RES(i);
end

end
function[phi,phiD,phiDD]=TestFuns(x)
% Test Functions
% phi(1)=(1+x^4)/2
% phi(2)=(1-x^4)
% phi(3)=(1-x^2)
phi(1)=(1+x^4)/2;phiD(1)=2*x^3;phiDD(1)=6*x^2;
phi(2)=1-x^4;phiD(2)=-4*x^3;phiDD(2)=-12*x^2;
%phi(2)=1-x^6;phiD(2)=-6*x^5;phiDD(2)=-30*x^4;
%phi(2)=cos(pi*x/2);phiD(2)=-(pi/2)*sin(pi*x/2);phiDD(2)=-(pi/2)^2*cos(pi*x/2);
phi(3)=1-x^2;phiD(3)=-2*x;phiDD(3)=-2;
end