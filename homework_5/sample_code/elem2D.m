function [ek,ef,node] = elem2D(XY,nel,k,dx,nodes,Case)
%Element stiffness and load matrices
%[ek][U1 U2 U3 U4]' = [f1 f2 f3 f4]'
%Case 1    -2nd Order Negative Case
%Case 2    -2nd Order Positive Case
%Case 3    -4th Order Negative Case
%Case 4    -4th Order Positive Case
%horiz_vert-denotes values that correspond to nodes above/below or
%           right/lift of node i,j in the stiffness matrix
%diagonal  -refers to values corresponding to nodes diagonal of node i,j in
%           the stiffness matrix
%center    -denotes the value corresponding to node i,j in stiffness matrix
%RHS       -denotes right hand side of fi


%Equations for each case are subject to depending on the differential
%equation. This code reflects +-(Uxx+Uyy)+k^2*U=k^2*x


switch Case 
    case 1 %2nd Order Negative Case
        
    %define molecule values
    horiz_vert = -1;
    center = (k^2*dx^2) + 4;
    RHS = k^2*dx^2;
    
    %elemental stiffness matrix
    ek = [center, 2*horiz_vert,0,2*horiz_vert;...
          2*horiz_vert,center,2*horiz_vert,0;...
          0,2*horiz_vert,center,2*horiz_vert;...
          2*horiz_vert,0,2*horiz_vert,center];
    
    %elemental load
    ef = [RHS*XY(nodes(nel,1),1);...
          RHS*XY(nodes(nel,2),1);...
          RHS*XY(nodes(nel,3),1);...
          RHS*XY(nodes(nel,4),1)];
  
    %define node numbers to be called later in assembly
    node = [nodes(nel,1);nodes(nel,2);nodes(nel,3);nodes(nel,4)];
 
    case 2  %2nd Order Positive Case
        
    %define molecule values
    horiz_vert = 1;
    center = (k^2*dx^2) - 4;
    RHS = k^2*dx^2;
    
    %elemental stiffness matrix
    ek = [center, 2*horiz_vert,0,2*horiz_vert;...
          2*horiz_vert,center,2*horiz_vert,0;...
          0,2*horiz_vert,center,2*horiz_vert;...
           2*horiz_vert,0,2*horiz_vert,center];
    
    %elemental load
    ef = [RHS*XY(nodes(nel,1),1);...
          RHS*XY(nodes(nel,2),1);...
          RHS*XY(nodes(nel,3),1);...
          RHS*XY(nodes(nel,4),1)];
  
%define node numbers to be called later in assembly
node = [nodes(nel,1);nodes(nel,2);nodes(nel,3);nodes(nel,4)];

    case 3 %4th Order Negative Case
        
    %define molecule values
    center = 40 + 8*k^2*dx^2;
    diagonal = -2;
    horiz_vert = k^2*dx^2 - 8;
    RHS = 12*k^2*dx^2;
    
    %elemental stiffness matrix
    ek = [center, 2*horiz_vert,4*diagonal,2*horiz_vert;...
          2*horiz_vert,center,2*horiz_vert,4*diagonal;...
          4*diagonal,2*horiz_vert,center,2*horiz_vert;...
          2*horiz_vert,4*diagonal,2*horiz_vert,center];
    
    %elemental load
    ef = [RHS*XY(nodes(nel,1),1);...
          RHS*XY(nodes(nel,2),1);...
          RHS*XY(nodes(nel,3),1);...
          RHS*XY(nodes(nel,4),1)];

    %define node numbers to be called later in assembly
    node = [nodes(nel,1);nodes(nel,2);nodes(nel,3);nodes(nel,4)];
    
    case 4 %4th Order Positive Case
        
    %define molecule values
    center = 8*k^2*dx^2 -40;
    diagonal = 2;
    horiz_vert =8 + k^2*dx^2;
    RHS = 12*k^2*dx^2;
    
    %elemental stiffness matrix
    ek = [center, 2*horiz_vert,4*diagonal,2*horiz_vert;...
          2*horiz_vert,center,2*horiz_vert,4*diagonal;...
          4*diagonal,2*horiz_vert,center,2*horiz_vert;...
          2*horiz_vert,4*diagonal,2*horiz_vert,center];
    
    %elemental load
    ef = [RHS*XY(nodes(nel,1),1);...
          RHS*XY(nodes(nel,2),1);...
          RHS*XY(nodes(nel,3),1);...
          RHS*XY(nodes(nel,4),1)];

    %define node numbers to be called later in assembly
    node = [nodes(nel,1);nodes(nel,2);nodes(nel,3);nodes(nel,4)];



end