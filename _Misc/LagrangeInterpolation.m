function f = LagrangeInterpolation(x,xNode,fNode)
    
    [Ndata,Nnode] = size(xNode);
    I       = 1:Nnode;
    
    
    BaryWeight  = zeros(Ndata,Nnode);
	Numerator   = zeros(size(x));
    
    for k = 1:Nnode
        Denominator = ones(size(x));
        for m = I(I~=k)
            Denominator = Denominator .* (xNode(:,k) - xNode(:,m));
        end
        BaryWeight(:,k) = 1./Denominator;
    end
    
    Denominator = zeros(size(x));
    for k = 1:Nnode
        InterpWeight = BaryWeight(:,k)./(x - xNode(:,k));
        Numerator    = Numerator   + InterpWeight .* fNode(:,k);
        Denominator  = Denominator + InterpWeight              ;
    end
    
    f = Numerator ./ Denominator;
    
end