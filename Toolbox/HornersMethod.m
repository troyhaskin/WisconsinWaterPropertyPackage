function Values = HornersMethod(x,p,mu)

    if(nargin < 3)
        mu = [0,1];
    end

    xstar  = (x - mu(1)) / mu(2)    ;
    Values = zeros(size(x)) + p(1)	; 
    N      = length(p)              ;

    for k = 2:N
        Values = p(k) + Values .* xstar;
    end

end