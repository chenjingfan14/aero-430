clear all; close all; clc
%#ok<*SPRIX>
%#ok<*CTCH>
%#ok<*UNRCH>
%#ok<*CLALL>
%#ok<*SAGROW>

%% 2nd-Order

meshOrder   = 3:4;
meshDx      = 1./(2.^meshOrder);

genFigSurface   = true;
genFigContour   = true;
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
        
        clear A ARow b u U
        
        dx = 1/(2^order);
        n = 1/dx;
        
        A = assemble(n, k, 2);
        
        b = sparse((n-1)^2, 1);
        
        for i = (n-1)*(n-2)+1:(n-1)^2
            
            xi = (i-(n-1)*(n-2))*dx;
            b(i) = sin(pi*xi);
            
        end
        
        u = A\b;
        
        u = reshape(u, n-1, n-1)';
        U = zeros(n+1, n+1);
        U(2:end-1, 2:end-1) = u;
        
        x = linspace(0, 1, n+1);
        y = linspace(0, 1, n+1);
        
        [X, Y] = meshgrid(x, y);
        
        U(end, :) = sin(pi.*x);
        
        titleString = strcat('2nd-Order CDS for the 2D Orthotropic Laplacian with k=', ...
            num2str(k), ' for dx=(1/2)^', num2str(order));
        figureString = strcat('order_2_k_', num2str(k), '_dx_order_', num2str(order));
        
        if genFigSurface
            figSurface(X, Y, U, titleString, figureString);
        end
        
        if genFigContour
            figContour(X, Y, U, titleString, figureString);
        end
        
        [dim, ~] = size(U);
        avg = (dim+1)/2;
        
        [~, Lprime] = lagrange(dx, avg, U, 2);
        
        uavgfdm(rowID, colID) = U(avg, avg);
        uavgexact(rowID, colID) = sin(pi/2)*sinh(k*pi/2)/sinh(k*pi);
        uyfdm(rowID, colID) = 1 / dx * ((1+k^2)*U(end, avg) - U(end-1, avg) - k^2/2*U(end, avg-1) - k^2/2*U(end, avg+1));
        uyexact(rowID, colID) = k*pi*coth(k*pi);
        uylag(rowID, colID) = double(subs(Lprime, 1));
        
        clear xl
        
    end
    
end

relError = abs(uyexact-uyfdm) ./ abs(uyexact) * 100;
relErrorAvg = abs(uavgexact-uavgfdm) ./ abs(uavgexact) * 100; %#ok<NASGU>
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
    figAvg(2, meshDx, relErrorAvg);
end

if genFigROC
    figRoc(2, meshDx, relError);
end

if genFigLag
    figLag(2, meshDx, relErrorLag);
end

if tableSave
    tableRender(2, meshDx, roc, rocLag, rocAvg, uyfdm, uyexact, relError);
end