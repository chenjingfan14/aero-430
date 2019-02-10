clear all; close all; clc

meshOrder   = 1:28;
meshDx      = (1/2).^meshOrder;

order = 2;

ufdm = sparse([]);
uexact = sparse([]);

index = 0;

for k = [1 2 5 10 20 50]
    
    index = index + 1;
    
    for dx = meshDx
        
        if order == 2
            
            kbar = 1/(pi*dx) * acosh(k^2+1-k^2*cos(pi*dx));
            
            ufdmcenter = sin(pi*1/2)*sinh(kbar*pi*1)/sinh(kbar*pi);
            ufdmright = sin(pi*(1/2+dx))*sinh(kbar*pi*1)/sinh(kbar*pi);
            ufdmleft = sin(pi*(1/2-dx))*sinh(kbar*pi*1)/sinh(kbar*pi);
            ufdmlower = sin(pi*1/2)*sinh(kbar*pi*(1-dx))/sinh(kbar*pi);
            
            uyfdm = 1 / dx * ((1+k^2)*ufdmcenter - ufdmlower - k^2/2*ufdmleft - k^2/2*ufdmright);
            uyexact = k*pi*coth(k*pi);
            
        elseif order == 4
            
            alpha = 1/12*(k^2+1);
            beta = 1-2/12*(k^2+1);
            gamma = k^2-2/12*(k^2+1);
            delta = (-2+4/12)*(k^2+1);

            kbar = 1/(pi*dx) * acosh((-2*gamma*cos(pi*dx)-delta)/(4*alpha*cos(pi*dx)+2*beta));
            
            n = 1/dx;
            
            alpha = -k^2/6;
            beta = 1 + k^2/3;
            C = sparse(gallery('tridiag', n+1, alpha, beta, alpha));
            C(1, 1) = 1;
            C(1, 2) = 0;
            C(end, end-1) = 0;
            C(end, end) = 1;

            d = sparse(n+1, 1);

            for i = 2:n
                ufdmcenter = sin(pi*i*dx)*sinh(kbar*pi*1)/sinh(kbar*pi);
                ufdmlower = sin(pi*i*dx)*sinh(kbar*pi*(1-dx))/sinh(kbar*pi);
                d(i) = ufdmcenter/dx - ufdmlower/dx + (k^2*pi^2*dx)/2*ufdmcenter + (k^4*pi^4*dx^3)/24*ufdmcenter;
            end

            uy = C\d;
            uyfdm = uy(n/2+1);
            uyexact = k*pi*coth(k*pi);
            
        end
        
        err(-log(dx)/log(2)) = abs((uyfdm - uyexact) / uyexact);
        
    end
    
    uxerr(index, :) = err;
    
end

figure
hold on
set(gcf, 'Position', [1 1 624 550])

for i = 1:index
    plot(-meshDx, uxerr(i, :), 'linewidth', 2);
end
plot([-10^0 -10^-9], [eps eps], 'k:', 'linewidth', 2)
legend('k=1', 'k=2', 'k=5', 'k=10', 'k=20', 'k=50', 'Machine Error')
xlabel('-\Delta x'); ylabel('\epsilon_r_e_l')
title('2-D Orthotropic Diffusion Equation  |  First Derivative Convergence')


set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
grid on; grid minor;
box on;