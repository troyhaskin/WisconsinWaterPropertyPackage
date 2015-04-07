function varargout = FilterList(PredicateOrMask,varargin)
    
    Nfiltered = length(varargin);
    Filtered  = cellfun( @(v) v(PredicateOrMask,:),varargin,'UniformOutput',false);
    
    
    if (nargout >= Nfiltered)
        varargout                  = cell(1,nargout);
        varargout(1:Nfiltered)     = Filtered       ;
        varargout{Nfiltered+1:end} = []             ;

    elseif (nargout == 1)
        varargout{1} = Filtered;

    end

end