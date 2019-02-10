function plot_mesh(XY,nodes)
% =========================================================================
% PLOTS THE MESH GENERATED IN THE PREVIOUS STEP, ELEMENT BY ELEMENT, SO THAT
% THE USER CAN MAKE SURE EVERYTHING WAS DONE CORRECTLY.  THE EDGES OF THE
% ELEMENT ARE PLOTTED USING STRAIGHT LINES, WHICH MAY NOT BE A TRUE
% REPRESENTATION OF THE ACTUAL PHYSICAL MESH, BUT IT IS A GOOD
% APPROXIMATION.
% =========================================================================
  
  global nelem
  
  % ORDER IN WHICH NODAL COORIDINATES ARE PULLED FROM THE OVERALL GLOBAL
  % MATRIX, XY
  order = [1 2 3 4];
  N = 4;
  hold on
  
  % PLOT EACH ELEMENT, AND LABEL THE ELEMENT
  % COLOR SIDE A DIFFERENT COLOR IF THEY ARE BOUNDARYS
  for nel=1:nelem
      for i=1:N
         X(i) = XY(nodes(nel,order(i)),1);
         Y(i) = XY(nodes(nel,order(i)),2);
      end
     
      % PLOT ONLY THE FIRST 9 POINTS OF X,Y (NODE 1 IS INCLUDED TWICE TO "CLOSE
      % THE LOOP"
      plot(X,Y)
      % THE ELEMENT NUMBER IS SHOWN AT THE LOCATION OF NODE 9 (ESSENTIALLY,
      % THE CENTER OF THE ELEMENT
      text(X(N),Y(N),num2str(nel));
  end
  hold off

end
