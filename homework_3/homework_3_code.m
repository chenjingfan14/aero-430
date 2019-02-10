clear all; close all; clc

%% Initial Conditions

plotGen     = false;
plotSave    = false;
tableSave   = false;
tableSave2  = false;
tableSave3  = false;

odeType = '1st Order US';
%odeType = '2nd Order CDS';
%odeType = '4th Order CDS';

odeString = '1st_Order_US';
%odeString = '2nd_Order_CDS';
%odeString = '4th_Order_CDS';

mesh.order = 1:18;
mesh.dx = 0.5.^mesh.order;

rowID = 0;

%% Boundary Value Problem Solution

for c = [1 2 5 10 20 50]
    
    rowID = rowID + 1;
    colID = 0;
    
    if plotGen
        
        fig1 = figure(1);
        xlabel('x');    ylabel('u(x)');
        grid on;        grid minor;
        box on;         hold on;
        set(gcf, 'Position', [1 1 624 550])
        
        if strcmpi(odeType, '2nd Order CDS')
            titleString = strcat('2nd Order CDS FDM for c=', num2str(c));
        elseif strcmpi(odeType, '4th Order CDS')
            titleString = strcat('4th Order CDS FDM for c=', num2str(c));
        elseif strcmpi(odeType, '1st Order US')
            titleString = strcat('1st Order US FDM for c=', num2str(c));
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
        
        if strcmpi(odeType, '2nd Order CDS')
            alpha = -1 - c*dx/2;
            beta = 2;
            gamma = -1 + c*dx/2;
        elseif strcmpi(odeType, '4th Order CDS')
            alpha = -1 - c*dx/2 - c^2*dx^2/12;
            beta = 2 + c^2*dx^2/6;
            gamma = -1 + c*dx/2 - c^2*dx^2/12;
        elseif strcmpi(odeType, '1st Order US')
            alpha = -1 - c*dx;
            beta = 2 + c*dx;
            gamma = -1;
        end
        
        A = gallery('tridiag', nx, alpha, beta, gamma);
        
        A(1, 1) = 1;    A(1, 2) = 0;        b(1) = 0;
        A(nx, nx) = 1;  A(nx, nx-1) = 0;    b(nx) = 1;
        b(2:nx-1) = 0;
        
        u = A\b;
        
        if plotGen && dx >= mesh.dx(8)
            
            figure(1)
            plot(x, u, 'linewidth', 1)
            
        end
        
        if strcmpi(odeType, '2nd Order CDS')
            dudx.fdm(rowID, colID) = (1 -c*dx/2)^-1 * (1 - u(end-1)) / dx;
        elseif strcmpi(odeType, '4th Order CDS')
            dudx.fdm(rowID, colID) = (1 -c*dx/2 + c^2*dx^2/6 - c^3*dx^3/24)^-1 * ...
                (1 - u(end-1)) / dx;
        elseif strcmpi(odeType, '1st Order US')
            dudx.fdm(rowID, colID) = (1 - u(end-1)) / dx;
        end
        
        ux.exact = (exp(c*x) - 1) / (exp(c) - 1);
        dudx.exact(rowID, colID) = c*exp(c) / (exp(c) - 1);
        
        if plotGen && dx >= mesh.dx(8)
            
            figure(2)
            plot(x(2:end-1), abs(ux.exact(2:end-1)'-u(2:end-1)) ./ ...
                abs(ux.exact(2:end-1)')*100, 'o-')
            
        end
        
        if tableSave2 && dx >= mesh.dx(5)
            
            colLabels = {'$x$', '$u(x)$', '$\bar{u}(x)$', '$\epsilon_{rel}$'};
            
            matrix2latex([x' u ux.exact' (abs(ux.exact'-u)./abs(ux.exact')*100)], ...
                strcat('pointwise_error_table_dx_', num2str(nx), '_', lower(odeString), '.tex'), ...
                'columnLabels', colLabels, 'alignment', 'c', 'format', '%1.2e')
            
        end
        
    end
    
    if plotGen
        
        figure(1)
        
        fplot(@(x) (exp(c*x)-1)/(exp(c) - 1), [0 1], '-k', 'linewidth', 1.5)
        
        legend('\Deltax = (1/2)^1', '\Deltax = (1/2)^2', '\Deltax = (1/2)^3', ...
            '\Deltax = (1/2)^4', '\Deltax = (1/2)^5', '\Deltax = (1/2)^6', ...
            '\Deltax = (1/2)^7', '\Deltax = (1/2)^8', 'Analytical Solution', ...
            'location', 'eastoutside')
        
        drawnow
        
        figure(2)
        
        legend('\Deltax = (1/2)^1', '\Deltax = (1/2)^2', '\Deltax = (1/2)^3', ...
            '\Deltax = (1/2)^4', '\Deltax = (1/2)^5', '\Deltax = (1/2)^6', ...
            '\Deltax = (1/2)^7', '\Deltax = (1/2)^8', 'location', 'eastoutside')
        set(gca, 'YScale', 'log')
        
        drawnow
        
        if plotSave
            
            figure(1)
            figureString = strcat('solution_', lower(odeString), '_k_', num2str(k));
            saveas(gcf, figureString, 'epsc')
            
            figure(2)
            figureString = strcat('pointwise_error_', lower(odeString), '_k_', num2str(c));
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
    
    for kID = 1:6
        loglog(-mesh.dx, relError(kID, :), '-o', 'linewidth', 1.25);
    end
    
    title(strcat('Convergence Plot for the ', odeType))
    legend('c=1', 'c=2', 'c=5', 'c=10', 'c=20', 'c=50')
    set(gca, 'XScale', 'log');  set(gca, 'YScale', 'log');
    drawnow
    
    if plotSave
        
        figureString = strcat('convergence_', lower(odeString));
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

colLabels = {'$\Delta x$', '$\beta(c=1)$', '$\beta(c=2)$', '$\beta(c=5)$', ...
    '$\beta(c=10)$', '$\beta(c=20)$', '$\beta(c=50)$'};

if tableSave
    
    matrix2latex([mesh.dx(1:end-1)' roc'], strcat('roc_table_', lower(odeString), '.tex'), ...
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
        strcat('qoi_1_table_', lower(odeString), '.tex'), 'columnLabels', colLabels1, ...
        'alignment', 'c', 'format', '%5.4f')
    matrix2latex([mesh.dx(1:end)' dudx.fdm(4:6,:)', dudx.exact(4:6, :)', relError(4:6, :)'], ...
        strcat('qoi_2_table_', lower(odeString), '.tex'), 'columnLabels', colLabels2, ...
        'alignment', 'c', 'format', '%5.4f')
    
end
