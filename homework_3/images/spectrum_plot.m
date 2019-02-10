clear all; close all; clc

%% Initial Conditions

ode.type = 'negative';
ode.order = 2;

%% Boundary Value Problem Solution

id = 1;

for c = linspace(10, 0, 200)
    
    dx = (1 / 2)^6;
    
    nx = 1 / dx + 1;
    x = linspace(0, 1, nx);
    
    A = zeros(nx, nx);
    b = zeros(nx, 1);
    
    alpha = - c/2*dx - 1;
    beta = 2;
    gamma = + c/2*dx - 1;
    
    for i = 1:nx
        if i == 1
            A(i, i) = 1;
            b(i) = 0;
        elseif i == nx
            A(i, i) = 1;
            b(i) = 1;
        else
            A(i, i-1) = alpha;
            A(i, i) = beta;
            A(i, i+1) = gamma;
            b(i) = 0;
        end
    end
    
    u(id, :) = A\b;
    
    id = id + 1;
end

X = x;
Y = linspace(10, 0, 200);
[X, Y] = meshgrid(X, Y)

figure
hold on
surf(X, Y, flipud(u), 'edgealpha', 0.5,'facealpha', 1)
%contour3(X, Y, flipud(u), 30, 'k-', 'linewidth', 2)
colormap spring
grid on;    grid minor
% hold on
% for w = 1:200
%     plot3(x, w/4*ones(length(x)), u(w, :), 'color', [0.3 0.8 0.1])
% end
% grid on; grid minor
% ylim([0 50]); zlim([0 1])
% %set(gca, 'YScale', 'log')
