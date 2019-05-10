function [xopt,fopt,X] = NelderMeade(fcn,xIni,varargin)
% [xopt,fopt] = NelderMeade(fcn,xIni,EPS)
%
% Nelder-Meade method to minimize function fcn. xIni can be a simplex or a 
% single point. The algorithm stops when the standard deviation of the 
% evaluation of the points is less than EPS (default eps). 
%
% Examples:
%
% fcn = @(x) (x(1)-3).^2/2 + (x(2)-2).^2;
% [xopt,fopt] = NelderMeade(fcn,[0 0.1;0.2 0.75; 1.1 0]')
%
% [xopt,fopt] = NelderMeade(fcn,[-1 -1]')

% M. Valenzuela
% March, 2019

if ~isempty(varargin)
   EPS = varargin{1};
else
   EPS = eps;
end

% Standard constants
alpha = 1;
gamma = 2;
rho = 0.5;
sigma = 0.5;

[n,m] = size(xIni); % number of dimensions

if m==1
   X = [xIni repmat(xIni,1,n)];
   for i=1:n
      if abs(xIni(i))>0.0001
         X(i,i+1) = X(i,i+1)*0.05;
      else
         X(i,i+1) = 0.000025;
      end
   end
else
   X = xIni;
end


f = zeros(1,n);  % evaluation of points

for i=1:n+1
   f(i) = fcn(X(:,i));
end

while 1
   [f,Ind] = sort(f);              % Sort
   if sqrt(std(f))<EPS             % stopping criteria
      break
   end
   X = X(:,Ind);
   x0 = mean(X(:,1:n),2);          % Centroid
   xr = x0 + alpha*(x0-X(:,n+1));  % Reflection point
   fxr = fcn(xr);
   if f(1)<=fxr && fxr<f(n+1)
      X(:,n+1) = xr;               % Reflection 1
      f(n+1) = fxr;
   elseif fxr<f(1)
      xe = x0 + gamma*(xr-x0);     % expansion point
      fxe = fcn(xe);
      if fxe < fxr
         X(:,n+1) = xe;            % Expansion
         f(n+1) = fxe;
      else
         X(:,n+1) = xr;            % Reflection 2
         f(n+1) = fxr;
      end
   else
      xc = x0 + rho*(X(:,n+1)-x0); % contraction point
      fxc = fcn(xc);
      if fxc<f(n+1)
         X(:,n+1) = xc;            % Contraction
         f(n+1) = fxc;
      else
         X(:,2:end) = X(:,1) + sigma*(X(:,2:end)-X(:,1)); % Shrink
         for i=2:n+1
            f(i) = fcn(X(:,i));
         end
      end
   end
end
xopt = X(:,1);
fopt = f(1);

end

