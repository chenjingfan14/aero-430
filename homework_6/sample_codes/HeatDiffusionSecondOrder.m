clc
clear

% exact solution
k=20;
del=(1/2)^(7);
    numpoints=1+1/del;
    numnodes=(numpoints-2)*(numpoints-2);
 coor=zeros(numnodes,2);
 xex=0:del:1;
 yex=0:del:1;
 temp=1;
 for c6=1:numpoints-2
     for c7=1:numpoints-2
         coor(temp,1)=del*c7;
        coor(temp,2)=del*c6;
        temp=temp+1;
    end
 end
coor;
lim=50;
uexact=zeros(numnodes,1);
for outer=1:numnodes
    xi=coor(outer,1);
    yi=coor(outer,2);
    sum=0;
    
    for c8=1:lim
        for c9=1:lim
            
            bnm=(k^2)*4*(xi)*((1-cos(c8*pi))/(c8*pi))*(1-cos(c9*pi))/(c9*pi);
             sum=sum+(bnm/(((c8)^2+(c9^2))*pi^2+k^2))*sin(c8*pi*xi)*sin(c9*pi*yi);
        end
    end
     uexact(outer,1)=sum;
    
end
 
 zex=zeros(numpoints,numpoints);
 tempnum=1;
    for c2=(numpoints-1):-1:2
        for c3=2:(numpoints-1)
        zex(c2,c3)=uexact(tempnum);
        tempnum=tempnum+1;
        end
    end
 
      zex=zex;  
      
      umid=zex((numpoints+1)/2,(numpoints+1)/2);
      uprimeexact=(zex((numpoints+1)/2,numpoints)-zex((numpoints+1)/2,numpoints-1))/del;
      
figure(1)
subplot(2,1,1)
      surf(xex,yex,zex)
      xlabel('X')
        ylabel('Y')
        zlabel('u')
    title('Exact Solution k=20')
    colorbar('location','eastoutside');
    
 

subplot(2,1,2)
contourf(xex,yex,zex)
xlabel('X')
        ylabel('Y')
         colorbar('location','eastoutside');
    
    umidexact = uexact(ceil(end/2), :);











%2nd order approximate  
k=20;
xlength = 1;
ylength = 1;

deltaX = [0.5^1,0.5^2,0.5^3,0.5^4,0.5^5,0.5^6,0.5^7];
deltaY=deltaX;
d=numel(deltaX);
error=zeros(1,d);
umidapproximatesecondorder=zeros(1,d);

for a=1:d
% 5x5 Matrix (CxD, All Side Edges of Matrix are Known Boundary Conditions
C=(xlength/deltaX(a))+1;
D=(ylength/deltaY(a))+1;

x=linspace(0,xlength,C);
y=linspace(0,ylength,D);

%Boundary Conditions 
BCBottom = 0;
BCLeft = 0;
BCRight = 0;
BCTop = 0;
BottomEdge = BCBottom*ones(1,D-2);
LeftEdge = BCLeft*ones(1,C-2);
RightEdge = BCRight*ones(1,C-2);
TopEdge = BCTop*ones(1,D-2);

N=C-2;
M=D-2;

%NumberOfUnknowns
NumberOfUnknowns = N*M;

%A Matrix in AT=B
A=zeros(N*M);

for i = 1:N*M
    A(i,i)=-(4/deltaX(a)^2)-k^2;
end 


for i=1:N-1
    for j=1:M
        A(i+(j-1)*N,i+(j-1)*N+1)=1/deltaX(a)^2;
        A(i+(j-1)*N+1,i+(j-1)*N)=1/deltaX(a)^2;
    end
end

for i=1:N
    for j=1:M-1
        A(i+(j-1)*N, i+j*N)=1/deltaY(a)^2;
        A(i+j*N, i+(j-1)*N)=1/deltaY(a)^2;
    end
end

%b Matrix
b =-(k^2)*ones(NumberOfUnknowns,1);
x=linspace(0,xlength,(xlength/deltaX(a))+1);
x(:,1)=[];
x(:,end)=[];
x=x';

r=kron(x,ones(N,1));

for i=1:NumberOfUnknowns
    b(i)=b(i)*r(i);
end 

%Solving for unknowns
T=A\b;

% Breaking up T matrix to plot in 3D
F= zeros(C,D);
%BoundaryCondition
F(end,:)= BCTop*ones(1,D);
s=1;
for j=1:N
    for i=1:N
        F(i+1,j+1)= T(s);
        s=s+1;
    end
end

x=linspace(0,xlength,(xlength/deltaX(a))+1);
[xx yy]=meshgrid(x,y);

figure(a+1)
subplot(2,1,1)
surf(xx, yy, F)
xlabel('Length Along Plate')
ylabel('Height Along Plate')
zlabel('u')
title(sprintf('Second Order Finite Difference for k=20,deltaX=deltaY=%d',deltaX(a)));
 colorbar('location','eastoutside');


subplot(2,1,2)
contourf(xx, yy, F)
xlabel('X')
        ylabel('Y')
         colorbar('location','eastoutside');

umidapproximatesecondorder(a) = T(ceil(end/2), :);
secondordererror(a)=abs(umidexact-umidapproximatesecondorder(a));




end

figure(d+2)
plot(-log(deltaX),log(secondordererror))
xlabel('-log(deltaX=deltaY)')
ylabel('log[uexact(1/2)-uapprox(1/2)]')
title('Second Order Error vs StepSize Log Plot')



umidapproximatesecondorder=umidapproximatesecondorder';
secondordererror=secondordererror';
Table1(:,1)=[deltaX];
Table1(:,2)=[umidapproximatesecondorder];
Table1(:,3)=[secondordererror];


%Richardson Second Order
secondorderurichard = zeros(1,d);
for i=2:d-1
    secondorderurichard(i)=(umidapproximatesecondorder(i)^2-umidapproximatesecondorder(i-1)*umidapproximatesecondorder(i+1))/(2*umidapproximatesecondorder(i)-umidapproximatesecondorder(i-1)-umidapproximatesecondorder(i+1));
end

secondorderrichardsonerror=zeros(1,d);
for i=2:d-1
    secondorderrichardsonerror(i)=abs(umidexact-secondorderurichard(i));
end


secondorderbetarichardson = zeros(1,d);
for i=2:d-1
secondorderbetarichardson(i)=(1/log(2))*log((secondorderurichard(i)-umidapproximatesecondorder(i-1))/(secondorderurichard(i)-umidapproximatesecondorder(i)));
end 

Table1(:,4)=[secondorderurichard'];     
Table1(:,5)=[secondorderrichardsonerror'];
Table1(:,6)=[secondorderbetarichardson'];
