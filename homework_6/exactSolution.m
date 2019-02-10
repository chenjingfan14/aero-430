function [ psum, nsum ] = exactSolution( k, xex, yex, nTerms )

psum = 0;
nsum = 0;

for n = 1:2:nTerms
    for m = 1:2:nTerms
        
        bnm = 4*k^2*xex*((1 - cos(n*pi))/(n*pi))*(1 - cos(m*pi))/(m*pi);        
        psum = psum + (bnm/(-((n)^2 + (m^2))*pi^2 + k^2))*sin(n*pi*xex)*sin(m*pi*yex);
        nsum = nsum + (bnm/(((n)^2 + (m^2))*pi^2 + k^2))*sin(n*pi*xex)*sin(m*pi*yex);
        
    end
end

end