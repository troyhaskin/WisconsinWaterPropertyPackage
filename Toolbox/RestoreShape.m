function varargout = RestoreShape(varargin)
   
    Nin       = length(varargin);
    Nout      = nargout         ;

    if(Nin ~= 2*Nout)
        error('Thermodynamics:RestoreShape:WrongInputOutputCount' ,...
              ['The number of output arguments must be exactly half ',...
              'of the input argument count.']);
    end
    
    for k = 1:2:Nout
        Array = varargin{ k };
        Size  = varargin{k+1};
        varargout{k} = reshape(Array,Size);
    end
    
end