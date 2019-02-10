clear all; close all; clc

%% Initial Conditions

ode.type = 'Positive';          % sign of the ODE ['positive' or 'negative']
ode.order = 4;                  % order of the FDM [2 or 4]

dudx.order = 2;                 % order of the extraction [1 or 2]

mesh.order = 1:6;               % range of dx values (i.e. 0.5^1 through 0.5^6 --> 1:6)
mesh.dx = 0.5.^mesh.order;

rowID = 0;

%% Boundary Value Problem Solution

for k = [1 2 5 10 20]
    
    figure
    xlabel('x');    ylabel('u(x)');
    grid on;        grid minor;
    box on;         hold on;
    
    if ode.order == 2
        titleString = strcat(ode.type, ' ODE with 2nd-Order FDM for k=', num2str(k));
    elseif ode.order == 4
        titleString = strcat(ode.type, ' ODE with 4th-Order FDM for k=', num2str(k));
    end
    
    title(titleString)
    
    rowID = rowID + 1;
    colID = 0;
    
    for dx = mesh.dx
        
        nx = 1 / dx + 1;
        x = linspace(0, 1, nx);
        
        A = zeros(nx, nx);
        b = zeros(nx, 1);
        
        if strcmpi(ode.type, 'positive') && ode.order == 2
            alpha = 1 / k^2 / dx^2;
            beta = -2 / k^2 / dx^2 + 1;
        elseif strcmpi(ode.type, 'positive') && ode.order == 4
            alpha = 1 / k^2 / dx^2 + 1/12;
            beta = -2 / k^2 / dx^2 + 10/12;
        elseif strcmpi(ode.type, 'negative') && ode.order == 2
            alpha = -1 / k^2 / dx^2;
            beta = 2 / k^2 / dx^2 + 1;
        elseif strcmpi(ode.type, 'negative') && ode.order == 4
            alpha = -1 / k^2 / dx^2 + 1/12;
            beta = 2 / k^2 / dx^2 + 10/12;
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
        
        u = A\b;
        
        plot(x, u, 'linewidth', 1)
        
        colID = colID + 1;
        
        if strcmpi(ode.type, 'positive') && dudx.order == 1
            dudx.fdm(rowID, colID) = - u(end-1) / dx;
        elseif strcmpi(ode.type, 'positive') && dudx.order == 2
            dudx.fdm(rowID, colID) = - u(end-1) / dx + dx * k^2 / 2;
        elseif strcmpi(ode.type, 'negative') && dudx.order == 1
            dudx.fdm(rowID, colID) = - u(end-1) / dx;
        elseif strcmpi(ode.type, 'negative') && dudx.order == 2
            dudx.fdm(rowID, colID) = - u(end-1) / dx - dx * k^2 / 2;
        end
        
        if strcmpi(ode.type, 'positive')
            dudx.exact(rowID, colID) = 1 - k * cos(k) / sin(k);
        elseif strcmpi(ode.type, 'negative')
            dudx.exact(rowID, colID) = 1 - k * cosh(k) / sinh(k);
        end
        
    end
    
    if strcmpi(ode.type, 'positive')
        fplot(@(x) x-sin(k*x)/sin(k), [0 1], '-k', 'linewidth', 1.5)
    elseif strcmpi(ode.type, 'negative')
        fplot(@(x) x-sinh(k*x)/sinh(k), [0 1], '-k', 'linewidth', 1.5)
    end
    
    legend('\Deltax = (1/2)^1', '\Deltax = (1/2)^2', '\Deltax = (1/2)^3', '\Deltax = (1/2)^4', ...
        '\Deltax = (1/2)^5', '\Deltax = (1/2)^6', 'Analytical Solution', ...
        'location', 'best')
    
    figureString = strcat(lower(ode.type), '_ode_order_', num2str(ode.order), '_k_', num2str(k));
    
    %saveas(gcf, figureString, 'epsc')
    
end

%% Convergence Analysis

figure
xlabel('-log_{10}(\Deltax)');   ylabel('log_{10}(\epsilon_{rel})');
grid on;                        grid minor;
box on;                         hold on;

logRelError = log10(abs(dudx.exact-dudx.fdm) ./ abs(dudx.exact) * 100);

for kID = 1:5
    plot(-log10(mesh.dx), logRelError(kID, :), '-o', 'linewidth', 1.25);
end

if ode.order == 2 && dudx.order == 1
    titleString = strcat(ode.type, ' ODE with 2nd-Order FDM - 1st-Order First Derivative Approximation');
elseif ode.order == 2 && dudx.order == 2
    titleString = strcat(ode.type, ' ODE with 2nd-Order FDM - 2nd-Order First Derivative Approximation');
elseif ode.order == 4 && dudx.order == 1
    titleString = strcat(ode.type, ' ODE with 4th-Order FDM - 1st-Order First Derivative Approximation');
elseif ode.order == 4 && dudx.order == 2
    titleString = strcat(ode.type, ' ODE with 4th-Order FDM - 2nd-Order First Derivative Approximation');
end

title(titleString)

legend('k=1', 'k=2', 'k=5', 'k=10', 'k=20')

figureString = strcat(lower(ode.type), '_ode_order_', num2str(ode.order), '_fd_order_', num2str(dudx.order));

%saveas(gcf, figureString, 'epsc')