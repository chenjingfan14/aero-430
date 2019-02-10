clear all; close all; clc

%% Initial Conditions

ode.type = 'negative';
ode.order = 2;

%% Boundary Value Problem Solution

id = 1;

for k = logspace(0, 1, 200)
    
    dx = (1 / 2)^5;
    
    nx = 1 / dx + 1;
    x = linspace(0, 1, nx);
    
    A = zeros(nx, nx);
    b = zeros(nx, 1);
    
    if strcmpi(ode.type, 'positive')
        alpha = 1 / k^2 / dx^2;
        beta = -2 / k^2 / dx^2 + 1;
    elseif strcmpi(ode.type, 'negative')
        alpha = -1 / k^2 / dx^2;
        beta = 2 / k^2 / dx^2 + 1;
    end
    
    for i = 1:nx
        if i == 1
            A(i, i) = 1;
            b(i) = 0;
        elseif i == nx
            A(i, i) = 1;
            b(i) = 0;
        else
            A(i, i-1) = alpha;
            A(i, i) = beta;
            A(i, i+1) = alpha;
            b(i) = x(i);
        end
    end
    
    u(id, :) = A\b;
    
    id = id + 1;
end

% if strcmpi(ode.type, 'positive')
%     fplot(@(x) x-sin(k*x)/sin(k), [0 1], '-k', 'linewidth', 1.5)
% elseif strcmpi(ode.type, 'negative')
%     fplot(@(x) x-sinh(k*x)/sinh(k), [0 1], '-k', 'linewidth', 1.5)
% end

figure
hold on
for w = 1:200
    plot3(x, w*ones(length(x)), u(w, :), 'color', [0.1 0.5 0.8])
end
grid on; grid minor
zlim([-10 10])
