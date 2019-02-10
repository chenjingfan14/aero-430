function [ ] = figLag( order, meshDx, relErrorLag )

figure
set(gcf, 'Position', [1 1 624 550])
xlabel('-log_{10}(\Deltax)');   ylabel('log_{10}(\epsilon_{rel})');
grid on;                        grid minor;
box on;                         hold on;

for kID = 1:6
    plot(-meshDx, relErrorLag(kID, :), '-o', 'linewidth', 1.25);
end

legend('1', '2', '5', '10', '20', '50', 'location', 'best')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

if order == 2
    title('2nd-Order CDS QOI Convergence using 2nd-Order Lagrangian Interpolation');
    saveas(gcf, 'order_2_u_y_lag', 'epsc');
elseif order == 4
    title('4th-Order CDS QOI Convergence using 2nd-Order/4th-Order Adaptive Lagrangian Interpolation');
    saveas(gcf, 'order_4_u_y_lag', 'epsc');
end


end

