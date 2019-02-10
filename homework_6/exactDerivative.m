function [ psum, nsum ] = exactDerivative( k, xex, yex, nTerms )

psum = 0;
nsum = 0;

for n = 1:2:nTerms
    for m = 1:2:nTerms
        
        bnm = 4*k^2*((1 - cos(n*pi))/(n*pi))*(1 - cos(m*pi))/(m*pi);
        psum = psum + sin(m*pi*yex)*(bnm/(-(n^2 + m^2)*pi^2 + k^2))*(sin(n*pi*xex)+n*pi*xex*cos(n*pi*xex));
        nsum = nsum + sin(m*pi*yex)*(bnm/((n^2 + m^2)*pi^2 + k^2))*(sin(n*pi*xex)+n*pi*xex*cos(n*pi*xex));
        
    end
end

end