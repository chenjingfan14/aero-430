function [ ] = figRoc( order, meshDx, relError )

figure
set(gcf, 'Position', [1 1 624 550])
xlabel('-log_{10}(\Deltax)');   ylabel('log_{10}(\epsilon_{rel})');
grid on;                        grid minor;
box on;                         hold on;

for kID = 1:6
    plot(-meshDx, relError(kID, :), '-o', 'linewidth', 1.25);
end

legend('1', '2', '5', '10', '20', '50', 'location', 'best')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

if order == 2
    title('2nd-Order CDS QOI Convergence using 2nd-Order FDM');
    saveas(gcf, 'order_2_u_y_fdm', 'epsc');
elseif order == 4
    title('4th-Order CDS QOI Convergence using 4th-Order FDM');
    saveas(gcf, 'order_4_u_y_fdm', 'epsc');
end


end