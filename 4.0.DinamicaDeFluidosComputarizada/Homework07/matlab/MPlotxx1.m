% Function of the vector field
function[X, Y, Ux, Uy] = MPlotxx1()
    [X,Y] = meshgrid(0.01 : 0.05 : 1.5, 0.01 : 0.05 : 1);
    Ux    = 2 * X .* Y;
    Uy    = -Y .^ 2;
    quiver(X, Y, Ux, Uy, 2, 'DisplayName', sprintf('gradient'));
    phi   = (Y .^ 3) / 3 - (X .^ 2) .* Y;
    contour(X, Y, phi, 25, 'DisplayName', sprintf('contours'));
end