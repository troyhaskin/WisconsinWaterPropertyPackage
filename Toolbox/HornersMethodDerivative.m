function varargout = HornersMethodDerivative(x,p,mu)

    if(nargin < 3)
        mu = [0,1];
    end

    xstar  = (x - mu(1)) / mu(2)    ;
    N      = length(p)              ;
    
    f    = zeros(size(x)) + p(1)    ; 
    dfdx = zeros(size(x)) + p(1)    ; 

    for k = 2:N-1
        f    = p(k) + f    .* xstar  ;
        dfdx = f    + dfdx .* xstar  ;
    end
    dfdx = dfdx/mu(2)         ; % Chain Rulle Correction arising from a non-zero scale in mu
    
    if (nargout == 1)
        varargout{1} = dfdx ;
    elseif(nargout == 2)
        
        f = p(N) + f .* xstar  ; % Complete the computation of f

        varargout{1} = f    ;
        varargout{2} = dfdx ;
    end
end