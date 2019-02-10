function [ ] = figContour( X, Y, U, titleString, figureString )

figure
contourf(X, Y, U)
colorbar
set(gcf, 'Position', [1 1 624 550])
title(titleString)
xlabel('x');    ylabel('y');    zlabel('u')
box on
%saveas(gcf, ['contour_' figureString], 'epsc')

end
