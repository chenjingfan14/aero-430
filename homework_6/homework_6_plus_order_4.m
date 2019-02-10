clear all; close all; clc
%#ok<*SPRIX>
%#ok<*CTCH>
%#ok<*UNRCH>
%#ok<*CLALL>
%#ok<*SAGROW>

%% 4th-Order

meshOrder   = 1:8;
meshDx      = 1./(2.^meshOrder);

kVals       = [1 2 5 10 20 50];

genFigSurface   = 0;
genFigContour   = 0;
genFigAvg       = 0;
tableSave       = 0;

rowID = 0;
colID = 0;

matNumTotal = numel(meshOrder)*numel(kVals);
matNum = 0;

for k = kVals
    
    rowID = rowID + 1;
    colID = 0;
    
    for order = meshOrder
        
        matNum = matNum + 1;
        
        clc
        fprintf('solving matrix %i/%i\n', matNum, matNumTotal)
        
        colID = colID + 1;
        
        clear A ARow b u U C d ux
        
        dx = 1/(2^order);
        n = 1/dx;
        
        A = assemblePlus(dx, n, k, 4);
        
        xi = k^2*dx^2*linspace(0+dx, 1-dx, n-1)';
        
        for i = 1:(n-1)
            if i == 1
                b = xi;
            else
                b = [b; xi];
            end
        end
        
        u = A\b;
        
        u = reshape(u, n-1, n-1)';
        U = zeros(n+1, n+1);
        U(2:end-1, 2:end-1) = u;
        
        x = linspace(0, 1, n+1);
        y = linspace(0, 1, n+1);
        
        [X, Y] = meshgrid(x, y);
        
        titleString = strcat('4th-Order CDS for the 2D Wave Equation with k=', ...
            num2str(k), ' for dx=(1/2)^', num2str(order));
        figureString = strcat('plus_order_4_k_', num2str(k), '_dx_order_', num2str(order));
        
        if genFigSurface
            figSurface(X, Y, U, titleString, figureString);
        end
        
        if genFigContour
            figContour(X, Y, U, titleString, figureString);
        end
        
        [dim, ~] = size(U);
        avg = (dim+1)/2;
        
        uavgfdm(rowID, colID) = U(avg, avg);
        [uavgexact(rowID, colID), ~] = exactSolution(k, 1/2, 1/2, 1000);
        
    end
    
end

relErrorAvg = abs(uavgexact-uavgfdm) ./ abs(uavgexact) * 100;
logRelErrorAvg = log10(relErrorAvg);

for kID = 1:numel(kVals)
    for rocID = 1:length(logRelErrorAvg) - 1
        rocAvg(kID, rocID) = (logRelErrorAvg(kID, rocID+1) - logRelErrorAvg(kID, rocID)) / -log10(2);
    end
end

if genFigAvg
    figAvg(4, meshDx, relErrorAvg);
end

if tableSave
    tableRender(4, meshDx, 0, rocAvg, 0, 0, 0, uavgfdm);
end