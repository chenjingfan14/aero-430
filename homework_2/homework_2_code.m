clear all; close all; clc

%% Initial Conditions

plotGen     = false;
plotSave    = false;
tableSave   = false;
tableSave2  = false;
tableSave3  = true;

odeType = 'Positive';
odeOrder = 4;

for dudxOrder = [4]
    
    mesh.order = 1:18;
    mesh.dx = 0.5.^mesh.order;
    
    rowID = 0;
    
    %% Boundary Value Problem Solution
    
    for k = [1 2 5 10 20 50]
        
        rowID = rowID + 1;
        colID = 0;
        
        if plotGen
            
            fig1 = figure(1);
            xlabel('x');    ylabel('u(x)');
            grid on;        grid minor;
            box on;         hold on;
            set(gcf, 'Position', [1 1 624 550])
            
            if odeOrder == 2
                titleString = strcat(odeType, ' ODE with 2nd-Order FDM for k=', num2str(k));
            elseif odeOrder == 4
                titleString = strcat(odeType, ' ODE with 4th-Order FDM for k=', num2str(k));
            end
            
            title(titleString)
            
            fig2 = figure(2);
            xlabel('x');    ylabel('\epsilon_{rel}');
            grid on;        grid minor;
            box on;         hold on;
            set(gcf, 'Position', [1 1 624 550])
            
            title(titleString)
            
        end
        
        for dx = mesh.dx
            
            nx = 1 / dx + 1;
            x = linspace(0, 1, nx);
            b = zeros(nx, 1);
            colID = colID + 1;
            
            if strcmpi(odeType, 'positive') && odeOrder == 2
                alpha = 1 / k^2 / dx^2;
                beta = -2 / k^2 / dx^2 + 1;
            elseif strcmpi(odeType, 'positive') && odeOrder == 4
                alpha = 1 / k^2 / dx^2 + 1/12;
                beta = -2 / k^2 / dx^2 + 10/12;
            elseif strcmpi(odeType, 'negative') && odeOrder == 2
                alpha = -1 / k^2 / dx^2;
                beta = 2 / k^2 / dx^2 + 1;
            elseif strcmpi(odeType, 'negative') && odeOrder == 4
                alpha = -1 / k^2 / dx^2 + 1/12;
                beta = 2 / k^2 / dx^2 + 10/12;
            end
            
            A = gallery('tridiag', nx, alpha, beta, alpha);
            
            A(1, 1) = 1;    A(1, 2) = 0;        b(1) = 0;
            A(nx, nx) = 1;  A(nx, nx-1) = 0;    b(nx) = 0;
            b(2:nx-1) = x(2:nx-1);
            
            u = A\b;
            
            if plotGen && dx >= mesh.dx(8)
                
                figure(1)
                plot(x, u, 'linewidth', 1)
                
            end
            
            if strcmpi(odeType, 'positive') && dudxOrder == 1
                dudx.fdm(rowID, colID) = - u(end-1) / dx;
            elseif strcmpi(odeType, 'positive') && dudxOrder == 2
                dudx.fdm(rowID, colID) = - u(end-1) / dx + dx * k^2 / 2;
            elseif strcmpi(odeType, 'positive') && dudxOrder == 4
                dudx.fdm(rowID, colID) = 1 / (1 - dx^2 * k^2 / 6) * (-u(end-1) / dx + dx * k^2 / 2 - ...
                    dx^2 * k^2 / 6 - dx^3 * k^4 / 24);
            elseif strcmpi(odeType, 'negative') && dudxOrder == 1
                dudx.fdm(rowID, colID) = - u(end-1) / dx;
            elseif strcmpi(odeType, 'negative') && dudxOrder == 2
                dudx.fdm(rowID, colID) = - u(end-1) / dx - dx * k^2 / 2;
            elseif strcmpi(odeType, 'negative') && dudxOrder == 4
                dudx.fdm(rowID, colID) = 1 / (1 + dx^2 * k^2 / 6) * (-u(end-1) / dx - dx * k^2 / 2 + ...
                    dx^2 * k^2 / 6 - dx^3 * k^4 / 24);
            end
            
            if strcmpi(odeType, 'positive')
                ux.exact = x - sin(k*x) ./ sin(k);
                dudx.exact(rowID, colID) = 1 - k * cos(k) / sin(k);
            elseif strcmpi(odeType, 'negative')
                dudx.exact(rowID, colID) = 1 - k * cosh(k) / sinh(k);
                ux.exact = x - sinh(k*x) ./ sinh(k);
            end
            
            if plotGen && dx >= mesh.dx(8)
                
                figure(2)
                plot(x(2:end-1), abs(ux.exact(2:end-1)'-u(2:end-1)) ./ ...
                    abs(ux.exact(2:end-1)')*100, 'o-')
                
            end
            
            if tableSave2 && dx >= mesh.dx(5)
                
                colLabels = {'$x$', '$u(x)$', '$\bar{u}(x)$', '$\epsilon_{rel}$'};
    
                matrix2latex([x' u ux.exact' (abs(ux.exact'-u)./abs(ux.exact')*100)], ...
                    strcat('dx_', num2str(nx), '_', lower(odeType), '_ode_', ...
                    num2str(odeOrder), '_fdm'), 'columnLabels', colLabels, ...
                    'alignment', 'c', 'format', '%1.2e')
                
            end
            
        end
        
        if plotGen
            
            figure(1)
            
            if strcmpi(odeType, 'positive')
                fplot(@(x) x-sin(k*x)/sin(k), [0 1], '-k', 'linewidth', 1.5)
            elseif strcmpi(odeType, 'negative')
                fplot(@(x) x-sinh(k*x)/sinh(k), [0 1], '-k', 'linewidth', 1.5)
            end
            
%             legend('\Deltax = (1/2)^1', '\Deltax = (1/2)^2', '\Deltax = (1/2)^3', ...
%                 '\Deltax = (1/2)^4', '\Deltax = (1/2)^5', '\Deltax = (1/2)^6', ...
%                 '\Deltax = (1/2)^7', '\Deltax = (1/2)^8', 'Analytical Solution', ...
%                 'location', 'eastoutside')
            
            drawnow
            
            figure(2)
            
%             legend('\Deltax = (1/2)^1', '\Deltax = (1/2)^2', '\Deltax = (1/2)^3', ...
%                 '\Deltax = (1/2)^4', '\Deltax = (1/2)^5', '\Deltax = (1/2)^6', ...
%                 '\Deltax = (1/2)^7', '\Deltax = (1/2)^8', 'location', 'eastoutside')
            ylim([10^-14 10^4])
            set(gca, 'YScale', 'log')
            
            drawnow
            
            if plotSave
                
                figure(1)
                figureString = strcat(lower(odeType), '_ode_order_', ...
                    num2str(odeOrder), '_k_', num2str(k));
                saveas(gcf, figureString, 'epsc')
                figure(2)
                figureString = strcat('error_', lower(odeType), '_ode_order_', ...
                    num2str(odeOrder), '_k_', num2str(k));
                saveas(gcf, figureString, 'epsc')
                
                close gcf;  close gcf
                
            end
            
        end
        
    end
    
    
    %% Convergence Analysis
    
    relError = abs(dudx.exact-dudx.fdm) ./ abs(dudx.exact) * 100;
    
    if plotGen
        
        figure
        xlabel('-log_{10}(\Deltax)');   ylabel('log_{10}(\epsilon_{rel})');
        grid on;                        grid minor;
        box on;                         hold on;
        ylim([10^-12 10^4])
        
        for kID = 1:5
            loglog(-mesh.dx, relError(kID, :), '-o', 'linewidth', 1.25);
        end
        
        if odeOrder == 2 && dudxOrder == 1
            titleString = strcat(odeType, ' ODE with 2nd-Order FDM - 1st-Order First Derivative Approximation');
        elseif odeOrder == 2 && dudxOrder == 2
            titleString = strcat(odeType, ' ODE with 2nd-Order FDM - 2nd-Order First Derivative Approximation');
        elseif odeOrder == 2 && dudxOrder == 4
            titleString = strcat(odeType, ' ODE with 2nd-Order FDM - 4th-Order First Derivative Approximation');
        elseif odeOrder == 4 && dudxOrder == 1
            titleString = strcat(odeType, ' ODE with 4th-Order FDM - 1st-Order First Derivative Approximation');
        elseif odeOrder == 4 && dudxOrder == 2
            titleString = strcat(odeType, ' ODE with 4th-Order FDM - 2nd-Order First Derivative Approximation');
        elseif odeOrder == 4 && dudxOrder == 4
            titleString = strcat(odeType, ' ODE with 4th-Order FDM - 4th-Order First Derivative Approximation');
        end
        
        title(titleString)
        legend('k=1', 'k=2', 'k=5', 'k=10', 'k=20')
        set(gca, 'XScale', 'log');  set(gca, 'YScale', 'log');
        drawnow
        
        if plotSave
            
            figureString = strcat(lower(odeType), '_ode_order_', num2str(odeOrder), ...
                '_fd_order_', num2str(dudxOrder));
            saveas(gcf, figureString, 'epsc')
            close gcf
            
        end
        
    end
    
    %% Rate of Convergence Analysis
    
    logRelError = log10(relError);
    
    for kID = 1:6
        
        for rocID = 1:length(logRelError) - 1
            roc(kID, rocID) = (logRelError(kID, rocID+1) - logRelError(kID, rocID)) / -log10(2);
        end
        
    end
    
    colLabels = {'$\Delta x$', '$\beta(k=1)$', '$\beta(k=2)$', '$\beta(k=5)$', ...
        '$\beta(k=10)$', '$\beta(k=20)$'};
    
    if tableSave
        
        matrix2latex([mesh.dx(1:17)' roc'], strcat(lower(odeType), '_ode_', ...
            num2str(odeOrder), '_fdm_', num2str(dudxOrder), '_dudx.tex'), ...
            'columnLabels', colLabels, 'alignment', 'c', 'format', '%5.4f')
        
    end
    
    colLabels1 = {'$\Delta x$', '$u''_{c=1}(1)$', '$u''_{c=2}(1)$', '$u''_{c=5}(1)$', ...
    '$\bar{u}''_{c=1}(1)$', '$\bar{u}''_{c=2}(1)$', '$\bar{u}''_{c=5}(1)$', ...
    '$\epsilon''_{rel,c=1}$', '$\epsilon''_{rel,c=2}$', '$\epsilon''_{rel,c=5}$'};
colLabels2 = {'$\Delta x$', '$u''_{c=10}(1)$', '$u''_{c=20}(1)$', '$u''_{c=50}(1)$', ...
    '$\bar{u}''_{c=10}(1)$', '$\bar{u}''_{c=20}(1)$', '$\bar{u}''_{c=50}(1)$', ...
    '$\epsilon''_{rel,c=10}$', '$\epsilon''_{rel,c=20}$', '$\epsilon''_{rel,c=50}$'};

if tableSave3
    
    matrix2latex([mesh.dx(1:end)' dudx.fdm(1:3,:)', dudx.exact(1:3, :)', relError(1:3, :)'], ...
        strcat('qoi_1_table_', lower(odeType),'_', num2str(odeOrder), '.tex'), 'columnLabels', colLabels1, ...
        'alignment', 'c', 'format', '%5.4f')
    matrix2latex([mesh.dx(1:end)' dudx.fdm(4:6,:)', dudx.exact(4:6, :)', relError(4:6, :)'], ...
        strcat('qoi_2_table_', lower(odeType),'_', num2str(odeOrder), '.tex'), 'columnLabels', colLabels2, ...
        'alignment', 'c', 'format', '%5.4f')
    
end
    
end