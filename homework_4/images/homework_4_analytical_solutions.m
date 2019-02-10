clc; clear all; close all

%%

figure
hold on

for k = [1 2 5 10 20 50]
%diff
fplot(@(x) x-sinh(k*x)/sinh(k), [0 1], 'linewidth', 1.5)
end

grid on;        grid minor
box on;
legend('k=1', 'k=2', 'k=5', 'k=10', 'k=20', 'k=50', 'location', 'best')
title('Analytical Solution to the Diffusion Equation')
xlabel('x');    ylabel('u(x)')

%%

figure
hold on

for k = [1 2 5 10 20 50]
%harm
fplot(@(x) x-sin(k*x)/sin(k), [0 1], 'linewidth', 1.5)
end

grid on;        grid minor
box on;
legend('k=1', 'k=2', 'k=5', 'k=10', 'k=20', 'k=50', 'location', 'best')
title('Analytical Solution to the Harmonic Wave Equation')
xlabel('x');    ylabel('u(x)')

%%

figure
hold on

for c = [1 2 5 10 20 50]
%conv
fplot(@(x) (exp(c*x)-1)/(exp(c) - 1), [0 1], 'linewidth', 1.5)
end

grid on;        grid minor
box on;
legend('c=1', 'c=2', 'c=5', 'c=10', 'c=20', 'c=50', 'location', 'best')
title('Analytical Solution to the Convection-Diffusion Equation')
xlabel('x');    ylabel('u(x)')