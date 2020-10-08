function[YOA,time,chi]=PelletPDE() 
    t0=0;
    tspan=linspace(0,2,100);
    [par]=Parameters();
    NP=par(1);
    bet=par(2);alf=par(3);m=par(5);C_inf=par(6);
    for i=1:NP
        y0(i)=0.0;
    end
    % Parameters for the orthogonal polynomials
    % ag=1,2,3 for planar, cylindrical and sphereical
    % be careful with al and be, if you do not know your problem
    al=1.0;be=0.5;ag=3.0;
    p=[al be ag];
    [time,Y]=ode45(@diffeqM,tspan,y0,[],p);
    plot(tspan,Y');
    YOA=Y';
    t=tspan;
    [xc,A,B,w,Q]=FJcbPlyx(al,be,ag,NP);
    [nrt,nct]=size(t);
    nrt=max(nrt,nct);
    Y=Y'; %may be not
    YT=Y;
    % The boundary condition can be applied to calcualte the
    %Concentration at the surface (i.e. dependent variable)
    % YT is the solution including the point at the surface
    for k=1:nrt
        sum3=0;
        for i=1:NP
            sum3=sum3+A(NP+1,i)*Y(i,k);
        end     
        YT(NP+1,k)=(bet*C_inf-sum3)/(bet+A(NP+1,NP+1));
    end
    % Concentration at the center, d are the coefficients of the 
    % polynomial for the solution
    d=inv(Q)*YT;
    % The solution will be extended to include the concentration
    % at the center
    % Y is the solution of the internal points
    % YT is the solution including the surface
    % YAO is the solution , including surface and center
    %chi is the dimensionless position , end includes center and surface
    chi=[0,xc']';
    YOA(1,:)=d(1,:); 
    YOA=[YOA;YT];
    plot(tspan,YOA);
    % To calcualte the average concentration
        for k=1:nrt
        sum5=0.0;
           for i=1:NP+1
               sum5=sum5+w(i)*YT(i,k);
           end
        sum5=sum5/sum(w);
        Yave(k)=sum5;
        end
        plot(tspan,YOA,tspan,Yave,'bo');
end
function [xc,A,B,w,Q] = FJcbPlyx( alpha,beta1,ag,N )
%Function to calculate the discretization coefficients using
%Orthogonal Collocation approach
% By: Dr. Jose Luis Lopez Salinas
%Updated version: 01-05-2015
% Mechanical Engineering Department (ITESM, Campus Monterrey)
%Warning: This is a general function so beta and ag are related
% beta=(ag-2)/2
%Calculates the coefficients of Jacobi Polynomials
% And also the roots of the polynomial of order N
% Integral (z^beta*(1-z)^alpha)*Ji*Jn*dz=0 , 0<z<1
% i=0,1,2,...,n-1
% Ji(z)=Jacobi Polynomial of order "i"
% r are the roots of the polynomial
% J(z)=a0+ap(1)*z+ap(2)*z^2+ap(3)*z^3+...+ap(N)*z^N
% J(z)=c(1)*z^N+c(2)*z^(N-1)+c(3)*z^(N-2)+...+c(N)*z+c(N+1)
% For the diffusion problem
% ag=1,2 or 3 depending of the geometry
% 1 is for planar, 2 is for cylindrical and 3 for spherical geometry
% alpha=1, and beta=(ag-2)/2
% Then the polynomial has the form:
% J(x^2)=a0+ap(1)*x^2+ap(2)*x^4+ap(3)*x^6+...+ap(N)*x^(2*N)
% Integral (1-x^2)*Pi(x^2)*Pn(x^2)*x^(ag-1)*dx, 0<x<1
% In this case the differential of volume is 
% dV=x^(ag-1)*dx
%=========================================================
% Differential equation of the form
% dC/dt - E Div (Grad(C))=R(C)
% dC/dt - E Lap (C)=R(C)
% Boundary condition of the form
% Grad(C)=-Bi(C-Co) at r/R=1
% Grad(C)=0 at r/R=0, this is already satisfied with symmetrical constrain
% Div = Divergence operator
% Grad = Gradient
% Lap=Laplacian = Lap()=Div(Grad())
% Sum=Summation
% =========================================================
% These matrices are used to discretize partial differential equations
% or just differential equations
% dc/dt - E Div (Grad(c))=R(c)
% dc(j)/dt - E Sum (B(j,i) c(i)) = [Lc^2/(E Co)]R(c(j)) {i=1,2,3...,N+1}
% B(j,i) = is an element of Matrix B
% A(j,i) = is an element of Matrix A
% w(j)=weights in the discretization
% Int (f*x^(ag-1)*dx)=Sum (w(i)*f(x(i)))
% {i=1,2,3...,N+1}, and 0<x<1
% Boundary condition
% -Grad(c)=Biot*(c-cinf) at r/R=1, this is at the external boundary r/R=1
% After discretization
% Sum ( A(N+1,i) c(i))= - Biot*(c(N+1)-cinf) {i=1,2,3...,N+1}
% Here Lc=Characteristic length
% E = Diffusion coeffieint
% Co = Reference molar concentration
% R(c)=Reaction rate in units of , kmol /(m^3-s)
% The rest of differnetial operators are dimensionless
% t=time*E/(Lc^2)
% x=r/R
% To calculate dependent variable at the center you can use
% The first term of the column vector callculated as:
% d=[Q^(-1)]*y
% So Q can be included in the output parameters of the function in the form
% function [c,ap0,ap,r,rx,xc,A,B,w,Q] = FJcbPly( alpha,beta,ag,N )
% =====================================================
a=alpha;b=beta1;
gamma0=1.0;
gamma1(1)=N*(N+1+a+b)/(b+1);
for i=2:N
    gamma1(i)=((N-i+1)*(N+i+a+b)/((i)*(i+b)))*gamma1(i-1);
end
ap0=gamma0*(-1)^(N);
for i=1:N
    ap(i)=((-1)^(N-i))*gamma1(i);
end
for i=1:N
    c(i)=ap(N+1-i);
end
c(N+1)=ap0;
if N==0
    ap=0;
end
r=roots(c);
rx=sqrt(r);
[nv,mv]=size(rx);
% Here the vectors are column vectors
% xc are collocation points
for i=1:N
    xc(i,1)=rx(N+1-i);
end
xc(N+1,1)=1.0;
x=xc;
for i=1:N+1
    for j=1:N+1
        Q(j,i)=x(j)^(2*i-2);
        C(j,i)=(2*i-2)*x(j)^(2*i-3);
        D(j,i)=(2*i-2)*(2*i+ag-4)*x(j)^(2*i-4);
    end
end
%A is the First Derivative Matrix
% B is Laplacian Matrix
% w the weights 
A=C*Q^(-1);
B=D*Q^(-1);

for i=1:N+1
    f(i)=1.0/(2*i-2+ag);
end
w=f*Q^(-1);
end
 % the differential equatuon
function[yprime]=diffeqM(t,y,p)
    [par]=Parameters();
    NE=par(1);
    for i=1:NE
        CP(i)=y(i);
    end
    al=p(1);be=p(2);ag=p(3);
    bet=par(2);alf=par(3);m=par(5);C_inf=par(6);
    [xc,A,B,w,Q]=FJcbPlyx(al,be,ag,NE);
    su1=0.0;
    for i=1:NE
        su1=su1+A(NE+1,i)*CP(i);
    end
    su1=su1-bet*C_inf;
    CP(NE+1)=-su1/(bet+A(NE+1,NE+1));
    for i=1:NE
        su2=0.0;
        for k=1:NE+1
         su2=su2+B(i,k)*CP(k);
        end
    zprime(i)=su2-alf*CP(i)^m;
end
   % Look out , this is a bug of the SciLab not of my code
    yprime=real(zprime);
    yprime=yprime';
end
% parameters of the differential equation
function[par]=Parameters()
    N=5;
    beta1=1.0;alpha=1;m=1;C_inf=1;Nx=4;
    par=[N,beta1,alpha,Nx,m,C_inf];
end