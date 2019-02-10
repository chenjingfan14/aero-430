clc; clear all; close all;
xi = 0; xf = 1; k =50; I = 10;
elements = i; Nodes = elements + 1;

%u value and dU(1)
x = xi:1/50:xf;
u = x-sinh(k*x)/sinh(k);
dU_true =  1-k*cosh(k*1)/sinh(k);
j = 1;

err = [];
RE_dU = [];
B= [];
err_RE = [];
j = 1;
ek = zeros(2,2)
for iterations = 1:I
    elements =2^iterations
 N = elements 
gk = zeros(N+1); gf = zeros(N+1,1);
del_x(j) = (xf - xi)/N;
dx = del_x(j);
X = linspace(xi,xf,N+1);
%Build elemental matrices
    for i = 1:N  
%element stiffness
    ek(1,1) = 2 + k^2*dx^2;
    ek(1,2) = -2;
    ek(2,1) = -2;
    ek(2,2) = 2 + k^2*dx^2;

%element load
    ef = [k^2*X(i)*dx^2;k^2*X(i+1)*dx^2];

%global stiffness
    gk(i,i) =  gk(i,i)+ ek(1,1); 
    gk(i+1,i) = gk(i+1,i) + ek(2,1);  
    gk(i,i+1) = gk(i,i+1) + ek(1,2);  
    gk(i+1,i+1) = gk(i+1,i+1) + ek(2,2);

%global loads
    gf(i) = gf(i) + ef(1); 
    gf(i+1) = gf(i+1) + ef(2);

    end

 %apply bc via penalty method
gk(1,1) = gk(1,1)+1e20; gk(N+1,N+1) = gk(N+1,N+1)+1e20; 

%SOLVE FOR DEFLECTION
j = iterations
    U = gk\gf;

    if j <2
    U1 = U;
    X1 = X
        elseif j <3
       U2 = U;
       X2 = X;     
        elseif j < 4
       U3 = U;
       X3 = X;
     
        elseif j < 5
       U4 = U;
       X4 = X;

        elseif j < 6
       U5 = U;
       X5 = X;
        elseif j < 7
       U6 = U;
       X6 = X;
     
        else
end
           
dU2(j) = (U(N+1)-U(N))/del_x(j)-(del_x(j)/2)*k^2;
err(j) = abs((dU_true-dU2(j))/dU_true*100);



j = j+1
end
for j =2:length(del_x)-1
RE_dU2(j)  = (dU2(j)^2-dU2(j-1)*dU2(j+1))/(2*dU2(j)-dU2(j-1)-dU2(j+1));
end
err_RE = zeros(length(del_x),1)
%err_RE -  error of richardson extrapolated value
for j =2:length(del_x)-1
err_RE(j) = abs(RE_dU2(j) - dU_true)/abs(dU_true);
end

%B- rate of convergence
for j =2:length(del_x)-1
B(j) = (1/log10(2))*log10((RE_dU2(j)-dU2(j-1))/(RE_dU2(j)-dU2(j)));
end





plot(x,u,X1,U1,X2,U2,X3,U3,X4,U4,X5,U5,X6,U6)
xlabel('x'); ylabel('y')
title('Modified FDM 2nd Order. K=50 (Negative)')
legend('exactsolution','\Delta x = 1/2^0','\Delta x = 1/2^1','\Delta x = 1/2^2','\Delta x = 1/2^3','\Delta x = 1/2^4','\Delta x = 1/2^5','\Delta x = 1/2^6')
figure
loglog(del_x,-err_RE,'o-')
title('Error for Modified FDM 2nd Order. K=50 (Negative)')
ylabel('log(err(dU(x=1)))'); xlabel('log(\Delta x)')