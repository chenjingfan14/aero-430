clc; clear all; close all;
k = 50;
%exact Solution
[Xt,Yt,Ut,dUt,midt,xx, slice_exact] = exact(k,1);
J = 1;
%define cartesian shape(in this case, a 1x1 square)
xi = 0; xf = 1; yi = 0; yf = 1; iterations = 5;
%number of segments in x&y direction
for j = 1:iterations
Nx = 2^j; Ny = Nx;
%Generate Mesh, &nodes
[XY,nodes,dx,dy] = mesh_square(xi,yi,xf,yf,Nx,Ny);

%number of elements and nodes 
nelem = length(nodes); nnodes = length(XY);

%define size of global stiffness and load matrices
gk = zeros(nnodes);
gf = zeros(nnodes,1);


for nel = 1:nelem
 %build element matricies
 %1-negative 2nd order; 2-positive 2nd order; 3-negative 4th order; 4-positive 4th order
[ek,ef,node] = elem2D(XY,nel,k,dx,nodes,1);

%compile matrices
[gk,gf] = assemble2D(ek,ef,gk,gf,node);
end

%Apply BC(All Edges are wquivalent to Zero)
for i = 1:length(XY)
        if mod(XY(i,1),1) < 0.0001
            gk(i,i) = gk(i,i) + 1e20;
            gf(i) = gf(i);
        elseif mod(XY(i,2),1) < 0.000001
            gk(i,i) = gk(i,i) + 1e20;
            gf(i) = gf(i);
        else
            gk(i,i) = gk(i,i);
            gf(i) = gf(i);
        end
end

%solve U
lin_U = gk\gf;
n = Nx;
%U is currently linearized. It must be modified to fit mesh.
Ufull = zeros(n+1,n+1);
for i = 1:n+1
    for j = 1:n+1
        Ufull(i,j) = lin_U((i-1)*(n+1)+j);
    end
end
x = []; x(1) = 0;
for i= 1:n
x(i+1) = x(i) + dx;
end
[X,Y] = meshgrid(x);

%1-d  view of x at y=0.5
for i=1:n+1
    slice(i)=Ufull(1+(n)/2,i);
end

%store values for different mesh sizes
if J < 2
    X1 = X; Y1 = Y; U1 = Ufull; slice1 = slice; x1=x;
elseif J < 3
    X2 = X; Y2 = Y; U2 = Ufull; slice2 = slice; x2=x;
elseif J < 4
    X3 = X; Y3 = Y; U3 = Ufull; slice3 = slice;  x3=x;
elseif J < 5
    X4 = X; Y4 = Y; U4 = Ufull; slice4 = slice; x4=x;
elseif J < 6
    X5 = X; Y5 = Y; U5 = Ufull; slice5 = slice;  x5=x;   
end

%store midpoint,dx,&error for each iteration
mid(J) = Ufull(1+(n)/2,1+(n)/2);
delx(J) = dx;
err(J) = abs(mid(J) - midt)/abs(midt);
J = J + 1;
end

%RE_mid - Richardson Extrapolation of the midpoint
for j =2:length(delx)-1
RE_mid(j)  = (mid(j)^2-mid(j-1)*mid(j+1))/(2*mid(j)-mid(j-1)-mid(j+1));
end

%err_RE -  error of richardson extrapolated value
for j =2:length(delx)-1
err_RE(j) = abs(RE_mid(j) - midt)/abs(midt);
end

%B- rate of convergence
for j =2:length(delx)-1
B(j) = (1/log10(2))*log10((RE_mid(j)-mid(j-1))/(RE_mid(j)-mid(j)));
end

%Plot contours
figure

subplot(3,2,[1,2])
suptitle('Negative Case,K=50 Exact & 2nd Order')
contourf(Xt,Yt,Ut);
xlabel 'x'; ylabel 'y'; title 'Exact solution'


subplot(3,2,3)
contourf(X2,Y2,U2)
xlabel 'x'; ylabel 'y'; title '\Delta x=1/4'

subplot(3,2,4)
contourf(X3,Y3,U3)
xlabel 'x'; ylabel 'y'; title '\Delta x=1/8'

subplot(3,2,5)
contourf(X4,Y4,U4)
xlabel 'x'; ylabel 'y'; title '\Delta x=1/16'

subplot(3,2,6)
contourf(X5,Y5,U5)
xlabel 'x'; ylabel 'y'; title '\Delta x=1/32'



%Plot slices

figure
plot(xx,slice_exact,x2,slice2,x3,slice3,x4,slice4,x5,slice5)
xlabel 'x'
ylabel 'U'
title 'Slice of U evauluated at y = 0.5 when K=50 (Negative Case)'
legend('exact','\Delta x=1/4','\Delta x=1/8','\Delta x=1/16','\Delta x=1/32')

%plot error of midpoint

figure
 loglog(-delx,err,'-o')
 title('error of U(0.5,0.5) when k=50 (Negative, 2nd Order)'); xlabel('log(\Delta X)'); ylabel('log(error(U(0.5,0.5))');
