function [gk,gf] = assemble2D(ek,ef,gk,gf,node)

N=4;
for i=1:N
    ig = node(i);
    
    % INSERT INTO GLOBAL LOAD VECTOR
    gf(ig) = gf(ig) + ef(i);
    
    for j=1:N
        jg = node(j);
        
        % INSERT INTO GLOBAL STIFFNESS MATRIX
        gk(ig,jg) = gk(ig,jg) + ek(i,j);
        
    end
end


end