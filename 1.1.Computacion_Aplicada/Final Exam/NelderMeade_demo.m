function [xopt,fopt,X] = NelderMeade_demo(fcn,limites,xIni)
% [xopt,fopt] = NelderMeade_demo(fcn,limites,xIni)
%
% Demostracion del metodo Nelder-Meade para minimizar una funcion fcn de
% dos variables. Se deben dar los limites para la grafica de contorno de
% fcn. Ini puede ser el simplex inicial o un punto. 
%
% Ejemplos:
%
% fcn = @(x) (x(1)-3).^2/2 + (x(2)-2).^2;
% [xopt,fopt] = NelderMeade_demo(fcn,[-1 5;-1 5],[-1 -1]')
%
% [xopt,fopt] = NelderMeade_demo(fcn,[-1 5;-1 5],[0 0.1;0.2 0.75; 1.1 0]')

% M. Valenzuela
% 0ct, 2018

EPS = 0.001;

p1 = linspace(limites(1,1),limites(1,2));
p2 = linspace(limites(2,1),limites(2,2));
[P1,P2] = meshgrid(p1,p2);
n = length(p1);
m = length(p2);
G = zeros(n,m);
for i=1:length(p1)
   for j=1:length(p2)
      G(i,j) = fcn([P1(i,j) P2(i,j)]);
   end
end
contour(P1,P2,G,20,'Color',[0.2 0.7 0.7],'LineStyle','-')
title('Ejemplo de Nelder-Meade')
grid on
xlabel('x_1')
ylabel('x_2')

% Constantes tipicas
alpha = 1;
gamma = 2;
rho = 0.5;
sigma = 0.5;

[n,m] = size(xIni); % numero de dimensiones

% Simplex inicial
if m==1
   % generado a partir de un punto
   X = [xIni repmat(xIni,1,n)];
   for i=1:n
      if abs(xIni(i))>eps 
         X(i,i+1) = X(i,i+1)*0.05;
      else
         X(i,i+1) = 0.000025;
      end
   end
else
   % dado por el usuario
   X = xIni;
end

f = zeros(1,n);  % evaluaciones de los puntos

hold on
plot([X(1,:) X(1,1)] ,[X(2,:) X(2,1)],'-or')
hold off
grid on
xlabel('x_1')
ylabel('x_2')
Xant = X;

for i=1:n+1
   f(i) = fcn(X(:,i));
end


while 1
   [f,Ind] = sort(f);              % Ordenamiento
   if sqrt(std(f))<EPS             % criterio de terminacion
      break
   end
   %    f
   X = X(:,Ind);
   pause
   x0 = mean(X(:,1:n),2);          % Centroide
   xr = x0 + alpha*(x0-X(:,n+1));  % punto de reflexion
   fxr = fcn(xr);
   if f(1)<=fxr && fxr<f(n+1)
      X(:,n+1) = xr;               % Reflexion 1
      f(n+1) = fxr;
      str = 'reflexion 1';
   elseif fxr<f(1)
      xe = x0 + gamma*(xr-x0);     % punto de expansion
      fxe = fcn(xe);
      if fxe < fxr
         X(:,n+1) = xe;            % Expansion
         f(n+1) = fxe;
         str = 'expansion';
      else
         X(:,n+1) = xr;            % Reflexion 2
         f(n+1) = fxr;
         str = 'reflexion 2';
      end
   else
      xc = x0 + rho*(X(:,n+1)-x0); % punto de contraccion
      fxc = fcn(xc);
      if fxc<f(n+1)
         X(:,n+1) = xc;            % Contraccion
         f(n+1) = fxc;
         str = 'contraccion';
      else
         X(:,2:end) = X(:,1) + sigma*(X(:,2:end)-X(:,1)); % Encogimiento
         for i=2:n+1
            f(i) = fcn(X(:,i));
         end
         str = 'encogimiento';
      end
   end
   contour(P1,P2,G,20,'Color',[0.2 0.7 0.7],'LineStyle','-')
   title(str)
   grid on
   xlabel('x_1')
   ylabel('x_2')
   hold on
   plot([X(1,:) X(1,1)] ,[X(2,:) X(2,1)],'.-r',...
      [Xant(1,:) Xant(1,1)] ,[Xant(2,:) Xant(2,1)],'.-b',...
      x0(1),x0(2),'ok',xr(1),xr(2),'xk')
   hold off
   Xant = X;
end

xopt = X(:,1);
fopt = f(1);


