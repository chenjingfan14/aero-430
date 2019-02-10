function [ L, Lprime ] = lagrange( dx, avg, U, lorder )

x0 = 1;
x1 = 1-dx;
x2 = 1-2*dx;
x3 = 1-3*dx;
x4 = 1-4*dx;
x5 = 1-5*dx;

syms xl

if lorder == 2
    lp0 = ((xl-x1)*(xl-x2))/((x0-x1)*(x0-x2));
    lp1 = ((xl-x0)*(xl-x2))/((x1-x0)*(x1-x2));
    lp2 = ((xl-x0)*(xl-x1))/((x2-x0)*(x2-x1));
    
    L = U(end, avg)*lp0 + U(end-1, avg)*lp1 + U(end-2, avg)*lp2;
end

if lorder == 4
    try
        lp0 = ((xl-x1)*(xl-x2)*(xl-x3)*(xl-x4))/((x0-x1)*(x0-x2)*(x0-x3)*(x0-x4));
        lp1 = ((xl-x0)*(xl-x2)*(xl-x3)*(xl-x4))/((x1-x0)*(x1-x2)*(x1-x3)*(x1-x4));
        lp2 = ((xl-x0)*(xl-x1)*(xl-x3)*(xl-x4))/((x2-x0)*(x2-x1)*(x2-x3)*(x2-x4));
        lp3 = ((xl-x0)*(xl-x1)*(xl-x2)*(xl-x4))/((x3-x0)*(x3-x1)*(x3-x2)*(x3-x4));
        lp4 = ((xl-x0)*(xl-x1)*(xl-x2)*(xl-x3))/((x4-x0)*(x4-x1)*(x4-x2)*(x4-x3));
        
        L = U(end, avg)*lp0 + U(end-1, avg)*lp1 + U(end-2, avg)*lp2 + U(end-3, avg)*lp3 + U(end-4, avg)*lp4;
    catch
        lp0 = ((xl-x1)*(xl-x2))/((x0-x1)*(x0-x2));
        lp1 = ((xl-x0)*(xl-x2))/((x1-x0)*(x1-x2));
        lp2 = ((xl-x0)*(xl-x1))/((x2-x0)*(x2-x1));
        
        L = U(end, avg)*lp0 + U(end-1, avg)*lp1 + U(end-2, avg)*lp2;
    end
end

Lprime = diff(L);

end

