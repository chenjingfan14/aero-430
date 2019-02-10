
clc; close all; clear all;

%odeString = 'diffusion_p_1';
%odeString = 'harmonic_p_1';
odeString = 'convection_p_1';

%odeType = 'Diffusion';
%odeType = 'Harmonic';
odeType = 'Convection';

meshOrder   = 1:18;
meshDx      = 0.5.^meshOrder;

rowID = 0;

plotGen = true;
plotSave = true;
tableSaveUx = false;
tableSaveRoc = false;
tableSaveError = false;

%% FEM for p = 1

for k = [1 2 5 10 20 50]
    
    rowID = rowID + 1;
    colID = 0;
    
    if plotGen
        
        fig1 = figure(1);
        xlabel('x');    ylabel('u(x)');
        grid on;        grid minor;
        box on;         hold on;
        set(gcf, 'Position', [1 1 624 550])
        
        if strcmpi(odeType, 'diffusion')
            titleString = strcat('p=1 Galerkin FEM for the Diffusion Equation for k=', num2str(k));
        elseif strcmpi(odeType, 'harmonic')
            titleString = strcat('p=1 Galerkin FEM for the Harmonic Equation for k=', num2str(k));
        elseif strcmpi(odeType, 'convection')
            titleString = strcat('p=1 Galerkin FEM for the Convection-Diffusion Equation for c=', num2str(k));
        end
        
        title(titleString)
        
        fig2 = figure(2);
        xlabel('x');    ylabel('\epsilon_{rel}');
        grid on;        grid minor;
        box on;         hold on;
        set(gcf, 'Position', [1 1 624 550])
        
        title(titleString)
        
    end
    
    for dx = meshDx
        
        nx = 1 / dx + 1;
        x = linspace(0, 1, nx);
        f = zeros(nx, 1);
        colID = colID + 1;
        
        for i = 2:nx-1
            if strcmpi('diffusion', odeType)
                f(i) = k^2*dx^2/6*(x(i-1)+4*x(i)+x(i+1));
            elseif strcmpi('harmonic', odeType)
                f(i) = k^2*dx^2/6*(x(i-1)+4*x(i)+x(i+1));
            elseif strcmpi('convection', odeType)
                f(i) = 0;
            end
        end
        
        if strcmpi('diffusion', odeType)
            alpha = -1 + k^2*dx^2/6;
            beta  = 2 + 4*k^2*dx^2/6;
            gamma = -1 + k^2*dx^2/6;
        elseif strcmpi('harmonic', odeType)
            alpha = 1 + k^2*dx^2/6;
            beta  = -2 + 4*k^2*dx^2/6;
            gamma = 1 + k^2*dx^2/6;
        elseif strcmpi('convection', odeType)
            alpha = -1 - k*dx/2;
            beta  = 2;
            gamma = -1 + k*dx/2;
        end
        
        A = gallery('tridiag', nx, alpha, beta, gamma);
        
        A(1, 1) = 1;    A(1, 2) = 0;
        A(nx, nx) = 1;  A(nx, nx-1) = 0;
        
        if strcmpi('diffusion', odeType)
            f(1) = 0;       f(nx) = 0;
        elseif strcmpi('harmonic', odeType)
            f(1) = 0;       f(nx) = 0;
        elseif strcmpi('convection', odeType)
            f(1) = 0;       f(nx) = 1;
        end
        
        u = A\f;
        
        if plotGen && dx >= meshDx(8)
            
            figure(1)
            plot(x, u, 'linewidth', 1)
            
        end
        
        if strcmpi(odeType, 'diffusion')
            dudx.fem(rowID, colID) = - u(end-1) / dx - dx * k^2 / 2;
            dudx.exact(rowID, colID) = 1 - k * cosh(k) / sinh(k);
            ux.exact = x - sinh(k*x) ./ sinh(k);
        elseif strcmpi(odeType, 'harmonic')
            dudx.fem(rowID, colID) = - u(end-1) / dx + dx * k^2 / 2;
            ux.exact = x - sin(k*x) ./ sin(k);
            dudx.exact(rowID, colID) = 1 - k * cos(k) / sin(k);
        elseif strcmpi(odeType, 'convection')
            dudx.fem(rowID, colID) = (1 -k*dx/2)^-1 * (1 - u(end-1)) / dx;
            ux.exact = (exp(k*x) - 1) / (exp(k) - 1);
            dudx.exact(rowID, colID) = k*exp(k) / (exp(k) - 1);
        end
        
        if plotGen && dx >= meshDx(8)
            
            figure(2)
            plot(x(2:end-1), abs(ux.exact(2:end-1)'-u(2:end-1)) ./ ...
                abs(ux.exact(2:end-1)')*100, 'o-')
            
        end
        
        if tableSaveUx && dx >= meshDx(5)
            
            colLabels = {'$x$', '$u(x)$', '$\bar{u}(x)$', '$\epsilon_{rel}$'};
            
            matrix2latex([x' u ux.exact' (abs(ux.exact'-u)./abs(ux.exact')*100)], ...
                strcat('pointwise_error_table_dx_', num2str(nx), '_', odeString, '.tex'), ...
                'columnLabels', colLabels, 'alignment', 'c', 'format', '%1.2e')
            
        end
        
    end
    
    if plotGen
        
        figure(1)
        
        if strcmpi(odeType, 'diffusion')
            fplot(@(x) x-sinh(k*x)/sinh(k), [0 1], '-k', 'linewidth', 1.5)
        elseif strcmpi(odeType, 'harmonic')
            fplot(@(x) x-sin(k*x)/sin(k), [0 1], '-k', 'linewidth', 1.5)
        elseif strcmpi(odeType, 'convection')
            fplot(@(x) (exp(k*x)-1)/(exp(k) - 1), [0 1], '-k', 'linewidth', 1.5)
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
        set(gca, 'YScale', 'log')
        
        drawnow
        
        if plotSave
            
            figure(1)
            figureString = strcat('solution_', lower(odeString), '_k_', num2str(k));
            saveas(gcf, figureString, 'epsc')
            
            figure(2)
            figureString = strcat('pointwise_error_', lower(odeString), '_k_', num2str(k));
            saveas(gcf, figureString, 'epsc')
            
            close gcf;  close gcf
            
        end
        
    end
    
end

%% Convergence Analysis

relError = abs(dudx.exact-dudx.fem) ./ abs(dudx.exact) * 100;

if plotGen
    
    figure
    xlabel('-log_{10}(\Deltax)');   ylabel('log_{10}(\epsilon_{rel})');
    grid on;                        grid minor;
    box on;                         hold on;
    ylim([10^-12 10^4])
    
    for kID = 1:6
        loglog(-meshDx, relError(kID, :), '-o', 'linewidth', 1.25);
    end
    
    if strcmpi(odeType, 'diffusion')
        titleString = strcat(odeType, ' Equation Convergence for p=1');
    elseif strcmpi(odeType, 'harmonic')
        titleString = strcat(odeType, ' Equation Convergence for p=1');
    elseif strcmpi(odeType, 'convection')
        titleString = strcat(odeType, ' Equation Convergence for p=1');
    end
    
    title(titleString)
    
    if ~strcmpi(odeType, 'convection')
        legend('k=1', 'k=2', 'k=5', 'k=10', 'k=20', 'k=50')
    else
        legend('c=1', 'c=2', 'c=5', 'c=10', 'c=20', 'c=50')
    end
    
    set(gca, 'XScale', 'log');  set(gca, 'YScale', 'log');
    drawnow
    
    if plotSave
        
        figureString = strcat('convergence_', odeString);
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

if ~strcmpi(odeType, 'convection')
    colLabels = {'$\Delta x$', '$\beta(k=1)$', '$\beta(k=2)$', '$\beta(k=5)$', ...
        '$\beta(k=10)$', '$\beta(k=20)$', '$\beta(k=50)$'};
else
    colLabels = {'$\Delta x$', '$\beta(c=1)$', '$\beta(c=2)$', '$\beta(c=5)$', ...
        '$\beta(c=10)$', '$\beta(c=20)$', '$\beta(c=50)$'};
end

if tableSaveRoc
    
    matrix2latex([meshDx(1:17)' roc'], strcat('roc_', odeString, '.tex'), ...
        'columnLabels', colLabels, 'alignment', 'c', 'format', '%5.4f')
    
end

if ~strcmpi(odeType, 'convection')
colLabels1 = {'$\Delta x$', '$u''_{k=1}(1)$', '$u''_{k=2}(1)$', '$u''_{k=5}(1)$', ...
    '$\bar{u}''_{k=1}(1)$', '$\bar{u}''_{k=2}(1)$', '$\bar{u}''_{k=5}(1)$', ...
    '$\epsilon''_{rel,k=1}$', '$\epsilon''_{rel,k=2}$', '$\epsilon''_{rel,k=5}$'};
colLabels2 = {'$\Delta x$', '$u''_{k=10}(1)$', '$u''_{k=20}(1)$', '$u''_{k=50}(1)$', ...
    '$\bar{u}''_{k=10}(1)$', '$\bar{u}''_{k=20}(1)$', '$\bar{u}''_{k=50}(1)$', ...
    '$\epsilon''_{rel,k=10}$', '$\epsilon''_{rel,k=20}$', '$\epsilon''_{rel,k=50}$'};
else
colLabels1 = {'$\Delta x$', '$u''_{c=1}(1)$', '$u''_{c=2}(1)$', '$u''_{c=5}(1)$', ...
    '$\bar{u}''_{c=1}(1)$', '$\bar{u}''_{c=2}(1)$', '$\bar{u}''_{c=5}(1)$', ...
    '$\epsilon''_{rel,c=1}$', '$\epsilon''_{rel,c=2}$', '$\epsilon''_{rel,c=5}$'};
colLabels2 = {'$\Delta x$', '$u''_{c=10}(1)$', '$u''_{c=20}(1)$', '$u''_{c=50}(1)$', ...
    '$\bar{u}''_{c=10}(1)$', '$\bar{u}''_{c=20}(1)$', '$\bar{u}''_{c=50}(1)$', ...
    '$\epsilon''_{rel,c=10}$', '$\epsilon''_{rel,c=20}$', '$\epsilon''_{rel,c=50}$'};
end

if tableSaveError
    
    matrix2latex([meshDx(1:end)' dudx.fem(1:3,:)', dudx.exact(1:3, :)', relError(1:3, :)'], ...
        strcat('qoi_1_table_', lower(odeString), '.tex'), 'columnLabels', colLabels1, ...
        'alignment', 'c', 'format', '%5.4f')
    matrix2latex([meshDx(1:end)' dudx.fem(4:6,:)', dudx.exact(4:6, :)', relError(4:6, :)'], ...
        strcat('qoi_2_table_', lower(odeString), '.tex'), 'columnLabels', colLabels2, ...
        'alignment', 'c', 'format', '%5.4f')
    
end