function [Xt,Yt,Ut,dUt,midt,xx, slice_exact] = exact(k,Case)
%Case 1 -Negative Exact Solution Derived via Galerkin Method
%Case 2 -Positive Exact Solution Derived via Galerkin Method
%
%
switch Case
    case 1 %Negative Case
%exact solution
iterations = 6;
for n = 1:iterations
    for m = 1: iterations
        anm(n,m) = 4*k^2*(cos(n*pi)*(cos(m*pi)-1))/((m*n*pi^2)*((n^2+m^2)*pi^2+k^2));
    end
end

xx = 0:.01:1;

for i = 1:length(xx)
    for j = 1:length(xx)
        x = xx(i);
        y = xx(j);
        
    U(i,j) = anm(1,1)*sin(1*pi*x)*sin(1*pi*y) + anm(2,1)*sin(2*pi*x)*sin(1*pi*y)...
            +anm(3,1)*sin(3*pi*x)*sin(1*pi*y) + anm(4,1)*sin(4*pi*x)*sin(1*pi*y)...
            +anm(5,1)*sin(5*pi*x)*sin(1*pi*y) + anm(6,1)*sin(6*pi*x)*sin(1*pi*y)...
            +anm(1,3)*sin(1*pi*x)*sin(3*pi*y) + anm(2,3)*sin(2*pi*x)*sin(3*pi*y)...
            +anm(3,3)*sin(3*pi*x)*sin(3*pi*y) + anm(4,3)*sin(4*pi*x)*sin(3*pi*y)...
            +anm(5,3)*sin(5*pi*x)*sin(3*pi*y) + anm(6,3)*sin(6*pi*x)*sin(3*pi*y)...
            +anm(1,5)*sin(1*pi*x)*sin(3*pi*y) + anm(2,5)*sin(2*pi*x)*sin(5*pi*y)...
            +anm(3,5)*sin(3*pi*x)*sin(5*pi*y) + anm(4,5)*sin(4*pi*x)*sin(5*pi*y)...
            +anm(5,5)*sin(5*pi*x)*sin(5*pi*y) + anm(6,5)*sin(6*pi*x)*sin(5*pi*y);
    end
end
[Xt,Yt] = meshgrid(xx);
Ut=U';
midt = Ut(51,51)


%1-d  view of x at y=0.5
for i=1:length(xx)
    slice_exact(i)=Ut((1+length(xx))/2,i);
end

%evaluate du/dx at x=1 and y = 0.5
x = 1; y = 0.5;
dUt =     1*pi*anm(1,1)*cos(1*pi*x)*sin(1*pi*y) + 2*pi*anm(2,1)*cos(2*pi*x)*sin(1*pi*y)...
            +3*pi*anm(3,1)*cos(3*pi*x)*sin(1*pi*y) + 4*pi*anm(4,1)*cos(4*pi*x)*sin(1*pi*y)...
            +5*pi*anm(5,1)*cos(5*pi*x)*sin(1*pi*y) + 6*pi*anm(6,1)*cos(6*pi*x)*sin(1*pi*y)...
            +1*pi*anm(1,3)*cos(1*pi*x)*sin(3*pi*y) + 2*pi*anm(2,3)*cos(2*pi*x)*sin(3*pi*y)...
            +3*pi*anm(3,3)*cos(3*pi*x)*sin(3*pi*y) + 4*pi*anm(4,3)*cos(4*pi*x)*sin(3*pi*y)...
            +5*pi*anm(5,3)*cos(5*pi*x)*sin(3*pi*y) + 6*pi*anm(6,3)*cos(6*pi*x)*sin(3*pi*y)...
            +1*pi*anm(1,5)*cos(1*pi*x)*sin(3*pi*y) + 2*pi*anm(2,5)*cos(2*pi*x)*sin(5*pi*y)...
            +3*pi*anm(3,5)*cos(3*pi*x)*sin(5*pi*y) + 4*pi*anm(4,5)*cos(4*pi*x)*sin(5*pi*y)...
            +5*pi*anm(5,5)*cos(5*pi*x)*sin(5*pi*y) + 6*pi*anm(6,5)*cos(6*pi*x)*sin(5*pi*y);
        
    case 2 % Negative Case        
iterations = 6;
for n = 1:iterations
    for m = 1: iterations
        anm(n,m) = 4*k^2*(cos(n*pi)*(cos(m*pi)-1))/((m*n*pi^2)*(k^2-(n^2+m^2)*pi^2));
    end
end

xx = 0:.01:1;

for i = 1:length(xx)
    for j = 1:length(xx)
        x = xx(i);
        y = xx(j);
        
    U(i,j) = anm(1,1)*sin(1*pi*x)*sin(1*pi*y) + anm(2,1)*sin(2*pi*x)*sin(1*pi*y)...
            +anm(3,1)*sin(3*pi*x)*sin(1*pi*y) + anm(4,1)*sin(4*pi*x)*sin(1*pi*y)...
            +anm(5,1)*sin(5*pi*x)*sin(1*pi*y) + anm(6,1)*sin(6*pi*x)*sin(1*pi*y)...
            +anm(1,3)*sin(1*pi*x)*sin(3*pi*y) + anm(2,3)*sin(2*pi*x)*sin(3*pi*y)...
            +anm(3,3)*sin(3*pi*x)*sin(3*pi*y) + anm(4,3)*sin(4*pi*x)*sin(3*pi*y)...
            +anm(5,3)*sin(5*pi*x)*sin(3*pi*y) + anm(6,3)*sin(6*pi*x)*sin(3*pi*y)...
            +anm(1,5)*sin(1*pi*x)*sin(3*pi*y) + anm(2,5)*sin(2*pi*x)*sin(5*pi*y)...
            +anm(3,5)*sin(3*pi*x)*sin(5*pi*y) + anm(4,5)*sin(4*pi*x)*sin(5*pi*y)...
            +anm(5,5)*sin(5*pi*x)*sin(5*pi*y) + anm(6,5)*sin(6*pi*x)*sin(5*pi*y);
    end
end
[Xt,Yt] = meshgrid(xx);
Ut=U';
midt = U(51,51);


%1-d  view of x at y=0.5
for i=1:length(xx)
    slice_exact(i)=U((1+length(xx))/2,i);
end

%evaluate du/dx at x=1 and y = 0.5
x = 1; y = 0.5;
dUt =     1*pi*anm(1,1)*cos(1*pi*x)*sin(1*pi*y) + 2*pi*anm(2,1)*cos(2*pi*x)*sin(1*pi*y)...
            +3*pi*anm(3,1)*cos(3*pi*x)*sin(1*pi*y) + 4*pi*anm(4,1)*cos(4*pi*x)*sin(1*pi*y)...
            +5*pi*anm(5,1)*cos(5*pi*x)*sin(1*pi*y) + 6*pi*anm(6,1)*cos(6*pi*x)*sin(1*pi*y)...
            +1*pi*anm(1,3)*cos(1*pi*x)*sin(3*pi*y) + 2*pi*anm(2,3)*cos(2*pi*x)*sin(3*pi*y)...
            +3*pi*anm(3,3)*cos(3*pi*x)*sin(3*pi*y) + 4*pi*anm(4,3)*cos(4*pi*x)*sin(3*pi*y)...
            +5*pi*anm(5,3)*cos(5*pi*x)*sin(3*pi*y) + 6*pi*anm(6,3)*cos(6*pi*x)*sin(3*pi*y)...
            +1*pi*anm(1,5)*cos(1*pi*x)*sin(3*pi*y) + 2*pi*anm(2,5)*cos(2*pi*x)*sin(5*pi*y)...
            +3*pi*anm(3,5)*cos(3*pi*x)*sin(5*pi*y) + 4*pi*anm(4,5)*cos(4*pi*x)*sin(5*pi*y)...
            +5*pi*anm(5,5)*cos(5*pi*x)*sin(5*pi*y) + 6*pi*anm(6,5)*cos(6*pi*x)*sin(5*pi*y);
 %Center point
 midt = Ut(51,51);
 
 %slice of U
 for i=1:length(xx)
    slice_exact(i)=Ut((1+length(xx))/2,i);
 end
        
end
        
