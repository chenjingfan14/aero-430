function [XY,nodes,dx,dy] = mesh_square(xi,yi,xf,yf,Nx,Ny)
% ============================================================
% SETS UP A MESH FOR A PIPE WITH FINS (ONLY 1/16 IS MESHED)
% INPUTS:   xi --   initial x position
%           xf --   final x position
%           yi --   initial position
%           yf --   final y position
%           Nx --   NUMBER OF INTERVALS in x direction
%           Ny--   NUMBER OF INTERVALS in y direction
%
% OUTPUT:   XY --   COORDINATES OF ALL NODES
%           nodes -- MATRIX CONTAINING GLOBAL NODE NUMBER FOR 
%                   EACH LOCAL NODE NUMBER (1 THRU 4)
%           BC --   BOUNDARY CONDITION FLAG FOR EACH ELEMENT FACE
%                   (1 FOR INNER RADIUS, 2 FOR OUTER RADIUS PLUS FIN)
%                   IF BC(e,i) IS 0, THEN FACE i OF ELEMENT e IS AN
%                   INTERNAL FACE
%           MAT --  MATERIAL NUMBER FOR EACH ELEMENT
%           dX  --  x increment size
%           dY  --  y increment size
% ============================================================

Xint = (xf - xi)/(Nx);
 Yint = (yf - yi)/(Ny);
dx = Xint; dy = Yint;
  % DETERMINE TOTAL NUMBER OF NODES AND TOTAL NUMBER OF ELEMENTS
  % ------------------------------------------------------------
 
  global nelem nnodes
  nelem = (Nx)^2;
  nnodes = (Nx+1)*(Ny+1);
  XY = zeros(nnodes,2);

  nodes = zeros(4,nelem);
    npt = 1;
  for ny=0:Ny
      for nx=0:Nx
          XY(npt,:) = [xi+nx*Xint,yi+ ny*Yint];
          npt = npt+1;
      end
  end
  
  
for  i = 1:nelem
  if i < Nx +1 
  nodes(1,i) = i;
  else if mod((i-1)^2,nelem) <0.001
   nodes(1,i) = nodes(2,i-1) +1;
      else 
   nodes(1,i) = nodes(2,i-1);
  end
  end
  nodes(2,i) = nodes(1,i) + 1;
  nodes(3,i) = nodes(1,i) + Nx + 2;
  nodes(4,i) = nodes(3,i) - 1;
 
end
nodes = nodes';
%confirms mesh
plot_mesh(XY,nodes)

  %Example Element
  
   %              face 3
   %            4-------3
   %            |       |
   %    face 4  |       |  face 2
   %            |       |
   %            1-------2
   %              face 1
