function varargout = Columnify(varargin)
    
    Nin       = length(varargin);
    Nout      = nargout         ;
    
    
    if    (  Nin == Nout) || (Nout == 0)
        varargout = ColumnifyWithoutSizes(varargin,Nin);
        
    elseif(2*Nin == Nout)
        varargout = ColumnifyWithSizes   (varargin,Nin);
        
    else
        error('Thermodynamics:Columnify:WrongInputOutputCount'        ,...
              ['The number of outputs must be equal to or double the ',...
               'number of input arguments.']);
    end
    
end

function ColumnedArrays = ColumnifyWithoutSizes(Arrays,N)
    ColumnedArrays = cell(1,N);
    for k = 1:N
        ColumnedArrays{k} = Arrays{k}(:);
    end
end

function ColumnedAndSizes = ColumnifyWithSizes(Arrays,N)
    ColumnedAndSizes = cell(1,2*N);
    
    m = 1;
    for k = 1:N
        Array                = Arrays{k}    ;
        ColumnedAndSizes{ m } = Array(:)    ;
        ColumnedAndSizes{m+1} = size(Array) ;
        m = m + 2;
    end
end

