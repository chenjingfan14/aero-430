syms A B C D E F G H I deltaX k

Piece1 = ((A-2*B+C)/deltaX^2);

Piece2 = ((1/deltaX^4)*F)-(2/deltaX^4)*A+(1/deltaX^4)*G-(2/deltaX^4)*D+(4/deltaX^4)*B-(2/deltaX^4)*E+(1/deltaX^4)*H-(2/deltaX^4)*C+(1/deltaX^4)*I;

uxx = Piece1-((deltaX^2)/12)*((Piece1*k^2)-Piece2);

Piece3 = ((D-2*B+E)/deltaX^2);

uyy = Piece3-((deltaX^2)/12)*((Piece3*k^2)-Piece2);

EntireLeftHandSide = (uxx+uyy)-B*k^2;

%Right Hand Side is -x_i*k^2

LeftHandSide=simplify(EntireLeftHandSide)
