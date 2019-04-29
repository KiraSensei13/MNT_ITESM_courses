%% ************************************************************************
% * AUTHOR(S) :
% *     Bruno González Soria          (A01169284)
% *     Antonio Osamu Katagiri Tanaka (A01212611)
% *
% * FILENAME :
% *     HW03.m
% *
% * DESCRIPTION :
% *     Computación Aplicada (Ene 19 Gpo 1)
% *     Homework on Simulated Annealing
% *
% * NOTES :
% *     
% *
% * START DATE :
% *     25 Apr 2019
% ************************************************************************

% warning('off')
% clc;
% close all;

%% ************************************************************************
% Problem 1:
% Using simulated annealing in the global optimization toolbox, solve the
% traveling salesman problem for 70 cities with coordinates generated in the
% following way:
% N = 70;
% rng(123);
% coordinates = rand(N,2);
% Upload to Blackboard a pdf file that contains the folowing:
%   A script of your solution
%   A plot of the best route found

N = 70;
rng(123);
coordinates = rand(N,2);
distances = TSPtable(coordinates);
objFcnTSP(distances);

route0 = randperm(N)';
TSPplot(route0,coordinates,'b',2)
    xlabel('x')
    ylabel('y')
    title(sprintf('TSP with %d cities',N))
options = optimoptions(@simulannealbnd,'DataType','custom',...
    'AnnealingFcn',@TSPinversion);

r0 = randperm(N);
TSPplot(r0,coordinates,'-')
    title(sprintf('initial random route, cost=%f',objFcnTSP(r0)))

r = simulannealbnd(@objFcnTSP,r0,[],[],options);
TSPplot(r,coordinates,'-')
    title(sprintf('cost=%f',objFcnTSP(r)))

function distances = TSPtable(coordinates)
[nCities,nx] = size(coordinates);
if nx ~= 2
    error('coordinates should be an (nCities) x 2 matrix')
end
distances = zeros(nCities);
for i=1:nCities
    for j=i: nCities
        distances(i,j) = ...
            sqrt( (coordinates(i,1)-coordinates(j,1))^2 + ...
            (coordinates(i,2)-coordinates(j,2))^2);
        distances(j,i) = distances(i,j);
    end
end
end

function f = objFcnTSP(varargin)
% f = objFcnTSP(distances)
% Loads distance matrix.
% f = objFcnTSP(route)
% Evaluates a route given a distances matrix.
persistent distances
if isempty(varargin)
    clear distances
else
    [n,m] = size(varargin{1});
    if n==m
        distances = varargin{1};
    else
        route = varargin{1};
        n = length(route);
        % Initialize with distance between last a first city
        f = distances(route(n), route(1));
        for i=2:n
            % Add distance from city i to city i-1
            f = f + distances(route(i-1), route(i));
        end
    end
end
end

function neighbor = TSPinversion(optimValues,varargin)
route = optimValues.x;
n = length(route);
m1 = floor(rand*n)+1;
m2 = mod(floor(rand*(n-1))+m1, n)+1;
n1 = min([m1 m2]);
n2 = max([m1 m2]);
neighbor = route;
neighbor(n1:n2) = route(n2:-1:n1);
end

function TSPplot(x,coordinates)
x = coordinates(:,1); % x coordinates
y = coordinates(:,2); % y coordinates
% Now plot these points and make sure you add 1 term to return to the starting point
plot(route,y)
end