% Ross Alexander
% 2.10.18

clc; close all; clear all;

meshOrder   = 1:15;
meshDx      = 0.5.^meshOrder;

rowID = 0;

plotGen         = true;
plotSave        = true;
tableSaveRoc    = true;

%% 2nd-Order FDM

for K = [1E-3 1E-2 1E-1 1E+0 1E+1 1E+2 1E+3]
    
    % Enumerate constants
    b = 1;
    k2 = 1;
    k1 = K*k2;
    
    rowID = rowID + 1;
    colID = 0;
    
    if plotGen
        
        xlabel('x');    ylabel('u(x)');
        grid on;        grid minor;
        box on;         hold on;
        set(gcf, 'Position', [1 1 624 550])
        
        titleString = strcat('2nd-Order FDM for the Bimaterial Heat Convection Equation for K=', ...
            num2str(K));
        title(titleString)
        
    end
    
    for dx = meshDx
        
        % Develop mesh-specific values for FDM calculations
        nxi = 1 / dx - 1;               % number of interior nodes
        nxt = 1 / dx + 1;               % number of total nodes (includes boundary nodes)
        xi = linspace(dx, 1-dx, nxi);   % x-values at interior nodes
        xt = linspace(0, 1, nxt);       % x-values at total nodes
        f = zeros(nxi, 1);              % interior load vector
        colID = colID + 1;
        
        nxiHalf = (nxi+1)/2;            % index for central interior node (x = 1/2) rel. to interior
        nxtHalf = (nxt+1)/2;            % index for central interior node (x = 1/2) rel. to total
        
        % Enumerate upper matrix and lower matrix constants
        alpha1 = -k1;                       alpha2 = -k2;           % subdiagonal values
        beta1  = b*dx^2+2*k1;               beta2  = b*dx^2+2*k2;   % diagonal values
        gamma1 = -k1;                       gamma2 = -k2;           % superdiagonal values
            
        if dx ~= 1/2
            
            % Construct upper matrix and lower matrix
            A1 = gallery('tridiag', nxiHalf, alpha1, beta1, gamma1);    % upper matrix for k1
            A2 = gallery('tridiag', nxiHalf, alpha2, beta2, gamma2);    % lower matrix for k2
            
            % Construct central vector for node at x = 1/2 for k1 and k2
            a = zeros(1, nxi);
            a(nxiHalf - 1) = -k1;
            a(nxiHalf)     = b*dx^2+k1+k2;
            a(nxiHalf + 1) = -k2;
            
            % Assemble entire interior matrix by combining upper, central, and lower matrices
            A = [A1(1:end-1, :) zeros(nxiHalf-1, nxiHalf-1);
                a; 
                zeros(nxiHalf-1, nxiHalf-1) A2(2:end, :)];
            
            % Enumerate load values at boundary nodes
            f(1) = 0;       f(nxi) = 100*k2;
            
        else
            
            % If dx = 1/2 the solution is trivial
            A = b*dx^2+k1+k2;
            f = k2*100;
            
        end
        
        % Solve for u and append values at boundary nodes
        ui = A\f;
        ut = [0; ui; 100];
        
        if plotGen && dx >= meshDx(8)
            plot(xt, ut, 'linewidth', 1)
        end
        
        %% Quantities of Interest
        
        % Develop analytical solution and corresponding coefficients
        a1 = sqrt(b/k1);
        a2 = sqrt(b/k2);
        
        C = [cosh(0) sinh(0) 0 0;
            0 0 cosh(a2) sinh(a2);
            cosh(a1/2) sinh(a1/2) -cosh(a2/2) -sinh(a2/2);
            k1*a1*sinh(a1/2) k1*a1*cosh(a1/2) -k2*a2*sinh(a2/2) -k2*a2*cosh(a2/2)];
        
        d = [0
            100;
            0;
            0];
        
        aSolCoeff = C\d;
        
        % Calculate midpoint temperature and heat flux from analytical and FDM solutions
        u12.exact(rowID, colID) = aSolCoeff(1)*cosh(a1/2) + aSolCoeff(2)*sinh(a1/2);
        u12.fdm(rowID, colID)   = ut(nxtHalf);
        Q12.exact(rowID, colID) = -k1*(aSolCoeff(1)*a1*sinh(a1/2) + aSolCoeff(2)*a1*cosh(a1/2));
        Q12.fdm(rowID, colID)   = -k1*(ut(nxtHalf)-ut(nxtHalf-1))/dx - b*dx/2*ut(nxtHalf);      % Q-
%         Q12.fdm(rowID, colID)   = -k2*(ut(nxtHalf+1)-ut(nxtHalf))/dx + b*dx/2*ut(nxtHalf);    % Q+
        
    end
    
    if plotGen
        
%         legend('\Deltax = (1/2)^1', '\Deltax = (1/2)^2', '\Deltax = (1/2)^3', ...
%             '\Deltax = (1/2)^4', '\Deltax = (1/2)^5', '\Deltax = (1/2)^6', ...
%             '\Deltax = (1/2)^7', '\Deltax = (1/2)^8', ...
%             'location', 'eastoutside')
        
        drawnow
        
        if plotSave

            figureString = strcat('solution_K_', num2str(rowID));
            saveas(gcf, figureString, 'epsc')
            
            close gcf
            
        end
        
    end
    
end

%% Convergence Plotting

% Calculate relative error
relErroru12 = abs(u12.exact-u12.fdm) ./ abs(u12.exact) * 100;
relErrorQ12 = abs(Q12.exact-Q12.fdm) ./ abs(Q12.exact) * 100;

if plotGen
    
    %%% Midpoint Temperature
    
    fig2 = figure(2);
    xlabel('-log_{10}(\Deltax)');   ylabel('log_{10}(\epsilon_{rel})');
    grid on;                        grid minor;
    box on;                         hold on;
    ylim([10^-12 10^4])
    set(gcf, 'Position', [1 1 624 550])
    
    % Plot log-log relative error in midpoint temperature v. -mesh size
    for KID = 1:7
        loglog(-meshDx, relErroru12(KID, :), '-o', 'linewidth', 1.25);
    end
    
    titleString = strcat('Convergence for Midpoint Temperature');
    title(titleString)
    
    legend('K=1E-3', 'K=1E-2', 'K=1E-1', 'K=1E+0', 'K=1E+1', 'K=1E+2', 'K=1E+3')
    
    set(gca, 'XScale', 'log');  set(gca, 'YScale', 'log');
    drawnow
    
    %%% Midpoint Heat Flux
    
    fig3 = figure(3);
    xlabel('-log_{10}(\Deltax)');   ylabel('log_{10}(\epsilon_{rel})');
    grid on;                        grid minor;
    box on;                         hold on;
    ylim([10^-12 10^4])
    set(gcf, 'Position', [1 1 624 550])
    
    % Plot log-log relative error in midpoint heat flux v. -mesh size
    for KID = 1:7
        loglog(-meshDx, relErrorQ12(KID, :), '-o', 'linewidth', 1.25);
    end
    
    titleString = strcat('Convergence for Midpoint Heat Flux');
    title(titleString)
    
    legend('K=1E-3', 'K=1E-2', 'K=1E-1', 'K=1E+0', 'K=1E+1', 'K=1E+2', 'K=1E+3')
    
    set(gca, 'XScale', 'log');  set(gca, 'YScale', 'log');
    drawnow
    
    if plotSave
        
        figure(2)
        saveas(gcf, 'convergence_u12', 'epsc')
        
        figure(3)
        saveas(gcf, 'convergence_Q12', 'epsc')
        
    end
    
end

%% Rate of Convergence Determination

logRelErroru12 = log10(relErroru12);

for KID = 1:7
    
    for rocID = 1:length(logRelErroru12) - 1
        rocu12(KID, rocID) = (logRelErroru12(KID, rocID+1) - logRelErroru12(KID, rocID)) / -log10(2);
    end
    
end

if tableSaveRoc

    colLabelsRoc = {'$\Delta x$', '$\beta(K=0.001)$', '$\beta(K=0.01)$', '$\beta(K=0.1)$', ...
        '$\beta(K=1)$', '$\beta(K=10)$', '$\beta(K=100)$', '$\beta(K=1000)$'};
    matrix2latex([meshDx(1:length(meshOrder)-1)' rocu12'], 'roc_u12.tex', ...
        'columnLabels', colLabelsRoc, 'alignment', 'c', 'format', '%5.4f')

end

logRelErrorQ12 = log10(relErrorQ12);

for KID = 1:7
    
    for rocID = 1:length(logRelErrorQ12) - 1
        rocQ12(KID, rocID) = (logRelErrorQ12(KID, rocID+1) - logRelErrorQ12(KID, rocID)) / -log10(2);
    end
    
end

if tableSaveRoc

    colLabelsRoc = {'$\Delta x$', '$\beta(K=0.001)$', '$\beta(K=0.01)$', '$\beta(K=0.1)$', ...
        '$\beta(K=1)$', '$\beta(K=10)$', '$\beta(K=100)$', '$\beta(K=1000)$'};
    matrix2latex([meshDx(1:length(meshOrder)-1)' rocQ12'], 'roc_Q12.tex', ...
        'columnLabels', colLabelsRoc, 'alignment', 'c', 'format', '%5.4f')

end