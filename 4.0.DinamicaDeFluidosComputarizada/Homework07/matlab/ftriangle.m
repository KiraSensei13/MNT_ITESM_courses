% https://math.stackexchange.com/questions/1344690/is-it-possible-to-find-the-vertices-of-an-equilateral-triangle-given-its-center
function[x, y] = ftriangle(p3)
% p3 are the parameters of the triangle
% x0, y0 is the location of the center and a the apotema
% send parameters xo, yo and a
    x0 = p3(1);
    y0 = p3(2);
    a  = p3(3)*2;
    n  = round(a * 100); % points per side
    
    top_vertex_x = x0;
    top_vertex_y = y0 + sqrt(3) / 3 * a;
    
    right_vertex_x = x0 + a / 2;
    right_vertex_y = y0 - sqrt(3) / 6 * a;
    
    left_vertex_x = x0 - a / 2;
    left_vertex_y = y0 - sqrt(3) / 6 * a;
    
    % line from top_vertex to right_vertex
    [xx, yy] = lineFunction( ...
        top_vertex_x, top_vertex_y, ...
        right_vertex_x, right_vertex_y, n);
    for i = 1 : n
        % fprintf('%i \n', i)
        x(i) = xx(i);
        y(i) = yy(i);
    end
    
    % line from right_vertex to left_vertex
    count = 1;
    [xx, yy] = lineFunction( ...
        right_vertex_x, right_vertex_y, ...
        left_vertex_x, left_vertex_y, n);
    for i = n + 1 : n * 2
        % fprintf('%i \n', i)
        x(i) = xx(count);
        y(i) = yy(count);
        count = count + 1;
    end
    
    % line from left_vertex to top_vertex
    count = 1;
    [xx, yy] = lineFunction( ...
        left_vertex_x, left_vertex_y, ...
        top_vertex_x, top_vertex_y, n);
    for i = n * 2 + 1 : n * 3
        % fprintf('%i \n', i)
        x(i) = xx(count);
        y(i) = yy(count);
        count = count + 1;
    end

    % plot(x, y, 'ro-')
end