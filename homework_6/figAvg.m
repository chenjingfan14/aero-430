function [ ] = figAvg( order, meshDx, relErrorAvg )

figure
set(gcf, 'Position', [1 1 624 550])
xlabel('-log_{10}(\Deltax)');   ylabel('log_{10}(\epsilon_{rel})');
grid on;                        grid minor;
box on;                         hold on;

for kID = 1:6
    plot(-meshDx, relErrorAvg(kID, :), '-o', 'linewidth', 1.25);
end

legend('1', '2', '5', '10', '20', '50', 'location', 'best')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

if order == 2
    title('2nd-Order CDS U{(1/2, 1/2)} Local Point Convergence');
    saveas(gcf, 'order_2_u_avg', 'epsc');
elseif order == 4
    title('4th-Order CDS U{(1/2, 1/2)} Local Point Convergence');
    saveas(gcf, 'order_4_u_avg', 'epsc');
end


end