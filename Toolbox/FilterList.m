function varargout = FilterList(Predicate,varargin)
    
    Nfiltered = length(varargin);
    Filtered  = cellfun( @(v) v(Predicate,:),varargin,'UniformOutput',false);
    
    if (nargout == 1)  || (nargout == 0)
        varargout{1} = Filtered;
    end
    if (nargout == Nfiltered)
        varargout = cell(1,Nfiltered);
        
        for k = 1:Nfiltered
            varargout{k} = Filtered{k};
        end
    end
end