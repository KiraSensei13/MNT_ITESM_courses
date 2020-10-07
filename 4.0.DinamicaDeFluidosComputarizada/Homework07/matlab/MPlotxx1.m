% Function of the vector field
function[X, Y, Ux, Uy] = MPlotxx1()
    [X,Y] = meshgrid(0.01 : 0.05 : 1.5, 0.01 : 0.05 : 1);
    Ux    = 2 * X .* Y;
    Uy    = -Y .^ 2;
    quiver(X, Y, Ux, Uy, 2);
    phi   = (Y .^ 3) / 3 - (X .^ 2) .* Y;
    contour(X, Y, phi);
    % [m,n]=size(X);
    % for i=1:m
    %   for j=1:n
    %   	phi(i,j)=Y(i,j)^3/3-X(i,j)^2*Y(i,j);
    %   end
    % end
    % contour(X,Y,phi);
end