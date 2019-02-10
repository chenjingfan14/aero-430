function [ ] = figSurface( X, Y, U, titleString, figureString )

figure
h = surf(X, Y, U);
h.EdgeAlpha = 0.5;
set(gcf, 'Position', [1 1 624 550])
title(titleString)
xlabel('x');    ylabel('y');    zlabel('u')
box on
grid on;         grid minor
%zlim([-1E-1 1])
%saveas(gcf, ['surface_' figureString], 'epsc')

end

