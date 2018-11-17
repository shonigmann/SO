function [s,rw] = randomWalk(n,seed)
    
    rng(seed);
    s = randn(n,1);
    rw = s;
    
    for i = 2:n
        
        rw(i) = rw(i)+rw(i-1);
        
    end

end