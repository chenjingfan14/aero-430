function [ A ] = assemble( n, k, order )

ac = sparse(zeros(n-1, n-1));
ai = sparse(zeros(n-1, n-1));
az = sparse(zeros(n-1, n-1));

if order == 2
    
    for i = 1:(n-1)
        for j = 1:(n-1)
            
            if i == j
                ac(i, j) = 2+2*k^2;
                ai(i, j) = -1;
            elseif abs(i - j) == 1
                ac(i, j) = -k^2;
                ai(i, j) = 0;
            end
            
        end
    end
    
elseif order == 4
    
    for i = 1:(n-1)
        for j = 1:(n-1)
            
            if i == j
                ac(i, j) = (-2+1/3)*(k^2+1);
                ai(i, j) = 1-k^2/6-1/6;
            elseif abs(i - j) == 1
                ac(i, j) = k^2-k^2/6-1/6;
                ai(i, j) = k^2/12+1/12;
            end
            
        end
    end
    
else
    error('Invalid order.')
end

for i = 1:(n-1)
    
    for j = 1:(n-1)
        
        if j == 1 && i == j
            ARow = ac;
        elseif j == 1 && abs(i - j) == 1
            ARow = ai;
        elseif j == 1 && abs(i - j) > 1
            ARow = az;
        elseif i == j
            ARow = [ARow ac];
        elseif abs(i - j) == 1
            ARow = [ARow ai];
        elseif abs(i - j) > 1
            ARow = [ARow az];
        end
        
    end
    
    if i == 1
        A = [ARow];
    else
        A = [A; ARow];
    end
    
end

end