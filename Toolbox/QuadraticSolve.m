function [PlusRoot,MinusRoot] = QuadraticSolve(a,b,c,mu)
    
    if (nargin < 4)
        mu = [0,1];
    end

    d = b./a;
    Discriminant = sqrt(d.^2 - 4*c./a  )  ;

    PlusRoot  = (-d + Discriminant)/2;
    MinusRoot = (-d - Discriminant)/2;

    PlusRoot  = mu(1) + mu(2)*PlusRoot  ;
    MinusRoot = mu(1) + mu(2)*MinusRoot ;

end