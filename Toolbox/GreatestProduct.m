function Greatest = GreatestProduct(varargin)
    
    Greatest = GreatestCore(@(x1,x2)IsGreater(x1,x2),varargin{:});
    
end

function TrueFalse = IsGreater(Argument1,Argument2)
    TrueFalse = prod(abs(Argument1)) > prod(abs(Argument2));
end