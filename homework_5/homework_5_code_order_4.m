clear all; close all; clc
%#ok<*SPRIX>
%#ok<*CTCH>
%#ok<*UNRCH>
%#ok<*CLALL>
%#ok<*SAGROW>

%% 4th-Order

meshOrder   = 1:8;
meshDx      = 1./(2.^meshOrder);

genFigSurface   = false;
genFigContour   = false;
genFigROC       = false;
genFigAvg       = false;
genFigLag       = false;
tableSave       = true;

rowID = 0;
colID = 0;

matNumTotal = length(meshOrder)*6;
matNum = 0;

for k = [1 2 5 10 20 50]
    
    rowID = rowID + 1;
    colID = 0;
    
    for order = meshOrder
        
        matNum = matNum + 1;
        
        clc
        fprintf('solving matrix %i/%i\n', matNum, matNumTotal)
        
        colID = colID + 1;
        
        clear A ARow b u U C d uy
        
        dx = 1/(2^order);
        n = 1/dx;
        
        A = assemble(n, k, 4);
        
        b = sparse((n-1)^2, 1);
        
        for i = (n-1)*(n-2)+1:(n-1)^2
            
            xl = (i-1-(n-1)*(n-2))*dx;
            xc = (i-(n-1)*(n-2))*dx;
            xr = (i+1-(n-1)*(n-2))*dx;
            b(i) = -(k^2/12+1/12)*(sin(pi*xl)+sin(pi*xr))-(1-k^2/6-1/6)*sin(pi*xc);
            
        end
        
        u = A\b;
        
        u = reshape(u, n-1, n-1)';
        U = zeros(n+1, n+1);
        U(2:end-1, 2:end-1) = u;
        
        x = linspace(0, 1, n+1);
        y = linspace(0, 1, n+1);
        
        [X, Y] = meshgrid(x, y);
        
        U(end, :) = sin(pi.*x);
        
        titleString = strcat('4th-Order CDS for the 2D Orthotropic Laplacian with k=', ...
            num2str(k), ' for dx=(1/2)^', num2str(order));
        figureString = strcat('order_4_k_', num2str(k), '_dx_order_', num2str(order));
        
        if genFigSurface
            figSurface(X, Y, U, titleString, figureString);
        end
        
        if genFigContour
            figContour(X, Y, U, titleString, figureString);
        end
        
        [dim, ~] = size(U);
        avg = (dim+1)/2;
        
        alpha = -k^2/6;
        beta = 1 + k^2/3;
        C = sparse(gallery('tridiag', n+1, alpha, beta, alpha));
        C(1, 1) = 1;
        C(1, 2) = 0;
        C(end, end-1) = 0;
        C(end, end) = 1;
        
        d = sparse(n+1, 1);
        
        for i = 2:n
            d(i) = U(end, i)/dx - U(end-1, i)/dx + (k^2*pi^2*dx)/2*U(end, i) + (k^4*pi^4*dx^3)/24*U(end, i);
        end
        
        uy = C\d;
        
        [~, Lprime] = lagrange(dx, avg, U, 4);
        
        uavgfdm(rowID, colID) = U(avg, avg);
        uavgexact(rowID, colID) = sin(pi/2)*sinh(k*pi/2)/sinh(k*pi);
        uyfdm(rowID, colID) = uy(avg);
        uyexact(rowID, colID) = k*pi*coth(k*pi);
        uylag(rowID, colID) = double(subs(Lprime, 1));
        
        clear xl
        
    end
    
end

relError = abs(uyexact-uyfdm) ./ abs(uyexact) * 100;
relErrorAvg = abs(uavgexact-uavgfdm) ./ abs(uavgexact) * 100;
relErrorLag = abs(uyexact-uylag) ./ abs(uyexact) * 100;
logRelError = log10(relError);
logRelErrorLag = log10(relErrorLag);
logRelErrorAvg = log10(relErrorAvg);

for kID = 1:6
    
    try
        for rocID = 1:length(logRelError) - 1
            roc(kID, rocID) = (logRelError(kID, rocID+1) - logRelError(kID, rocID)) / -log10(2);
            rocLag(kID, rocID) = (logRelErrorLag(kID, rocID+1) - logRelErrorLag(kID, rocID)) / -log10(2);
            rocAvg(kID, rocID) = (logRelErrorAvg(kID, rocID+1) - logRelErrorAvg(kID, rocID)) / -log10(2);
        end
    catch
        for rocID = 1:length(logRelError) - 2
            roc(kID, rocID) = (logRelError(kID, rocID+1) - logRelError(kID, rocID)) / -log10(2);
            rocLag(kID, rocID) = (logRelErrorLag(kID, rocID+1) - logRelErrorLag(kID, rocID)) / -log10(2);
            rocAvg(kID, rocID) = (logRelErrorAvg(kID, rocID+1) - logRelErrorAvg(kID, rocID)) / -log10(2);
        end
    end
    
end

if genFigAvg
    figAvg(4, meshDx, relErrorAvg);
end

if genFigROC
    figAvg(4, meshDx, relError);
end

if genFigLag
    figLag(4, meshDx, relErrorLag);
end

if tableSave
    tableRender(4, meshDx, roc, rocLag, rocAvg, uyfdm, uyexact, relError);
end