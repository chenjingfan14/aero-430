clear all; close all; clc

set(0,'DefaultFigureWindowStyle','docked')

figure
cmap = colormap(hot);
cIndex = 1;

for K = [1E-3 1E-2 1E-1 1E0 1E1 1E2 1E3]
    
    b = 16;
    k2 = 1;
    k1 = K*k2;
    Ta = 0;
    A = 1;
    
    a1 = sqrt(b/k1);
    a2 = sqrt(b/k2);
    
    C = [cosh(0) sinh(0) 0 0;
        0 0 cosh(a2) sinh(a2);
        cosh(a1/2) sinh(a1/2) -cosh(a2/2) -sinh(a2/2);
        k1*a1*sinh(a1/2) k1*a1*cosh(a1/2) -k2*a2*sinh(a2/2) -k2*a2*cosh(a2/2)];
    
    d = [-Ta
        100-Ta;
        0;
        0];
    
    aSolCoeff = C\d;
    
    
    hold on;    box on;
    grid on;    grid minor;
    x1 = linspace(0.0, 0.5, 100);
    x2 = linspace(0.5, 1.0, 100);
    s1 = aSolCoeff(1)*cosh(a1*x1) + aSolCoeff(2)*sinh(a1*x1) + Ta;
    s2 = aSolCoeff(3)*cosh(a2*x2) + aSolCoeff(4)*sinh(a2*x2) + Ta;
    plot([x1 x2], [s1 s2],  'linewidth', 2, 'color', cmap(cIndex, :));
    xlim([0 1]);    ylim([0 100]);
    title('Analytical Temperature Profile')
    xlabel('Distance');  ylabel('Temperature')
    
    cIndex = cIndex + 9;
    
end

legend('K = 1E-3', 'K = 1E-2', 'K = 1E-1', 'K = 1E+0', 'K = 1E+1', 'K = 1E+2', 'K = 1E+3', ...
    'location', 'best') 
saveas(gcf, '430_hw1_analytical_temperature', 'epsc')

figure
cmap = colormap(hot);
cIndex = 1;

for K = [1E-3 1E-2 1E-1 1E0 1E1 1E2 1E3]
    
    b = 1;
    k2 = 1;
    k1 = K;
    Ta = 0;
    A = 1;
    
    a1 = sqrt(b/k1);
    a2 = sqrt(b/k2);
    
    C = [cosh(0) sinh(0) 0 0;
        0 0 cosh(a2) sinh(a2);
        cosh(a1/2) sinh(a1/2) -cosh(a2/2) -sinh(a2/2);
        k1*a1*sinh(a1/2) k1*a1*cosh(a1/2) -k2*a2*sinh(a2/2) -k2*a2*cosh(a2/2)];
    
    d = [-Ta
        100-Ta;
        0;
        0];
    
    aSolCoeff = C\d;
    
    hold on;    box on;
    grid on;    grid minor;
    x1 = linspace(0.0, 0.5, 100);
    x2 = linspace(0.5, 1.0, 100);
    s1 = -k1*A*(aSolCoeff(1)*a1*sinh(a1*x1) + aSolCoeff(2)*a1*cosh(a1*x1));
    s2 = -k2*A*(aSolCoeff(3)*a2*sinh(a2*x2) + aSolCoeff(4)*a2*cosh(a2*x2));
    plot([x1 x2], [s1 s2],  'linewidth', 2, 'color', cmap(cIndex, :));
    xlim([0 1]);    ylim([-inf 0]); 
    title('Analytical Heat Flux Profile')
    xlabel('Distance');  ylabel('Heat Flux')
    
    cIndex = cIndex + 9;
    
end

legend('K = 1E-3', 'K = 1E-2', 'K = 1E-1', 'K = 1E+0', 'K = 1E+1', 'K = 1E+2', 'K = 1E+3', ...
    'location', 'best') 
saveas(gcf, '430_hw1_analytical_heat_flux', 'epsc')