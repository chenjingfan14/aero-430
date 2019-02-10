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
        
        clear ufdm
        
        if order == 2
            
            kbar = 1 / dx * acosh(1 + k^2*dx^2/2);
            
%             x = 0:dx:1;
%             ufdm(-log(dx)/log(2), 1:1/dx+1) = x - sinh(kbar*x)/sinh(kbar);
%             uxfdm = -ufdm(-log(dx)/log(2), end-1) / dx - dx * k^2 / 2;
            
            ufdm = 1-dx - sinh(kbar*(1-dx))/sinh(kbar);
            uxfdm = -ufdm / dx - dx * k^2 / 2;
            
            uxexact = 1 - k*coth(k);
            
        elseif order == 4
            
            kbar = 1 / dx * acosh((2+10*k^2*dx^2/12) / (-2*(-1+k^2*dx^2/12)));
            
%             x = 0:dx:1;
%             ufdm(-log(dx)/log(2), 1:1/dx+1) = x - sinh(kbar*x)/sinh(kbar);
%             uxfdm = 1 / (1 + dx^2 * k^2 / 6) * (-ufdm(-log(dx)/log(2), end-1) / dx - dx * k^2 / 2 + ...
%                     dx^2 * k^2 / 6 - dx^3 * k^4 / 24);
            
            ufdm = 1-dx - sinh(kbar*(1-dx))/sinh(kbar);
            uxfdm = 1 / (1 + dx^2 * k^2 / 6) * (-ufdm / dx - dx * k^2 / 2 + ...
                    dx^2 * k^2 / 6 - dx^3 * k^4 / 24);        
                
            uxexact = 1 - k*coth(k);
            
        end
        
        err(-log(dx)/log(2)) = abs((uxfdm - uxexact) / uxexact);
        
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
title('1-D Diffusion Equation  |  First Derivative Convergence')

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
grid on; grid minor;
box on;