function [gK, gF] = BC(gk,gf,XY)
%apply BC via penalty method
for i = 1:length(XY)
        if mod(XY(i,1),1) < 0.00001
            gK(i,i) = gk(i,i) + 1e20;
            gF(i) = gf(i);
        else
            gK(i,i) = gk(i,i);
            gF(i) = gf(i);
        end
end

end
