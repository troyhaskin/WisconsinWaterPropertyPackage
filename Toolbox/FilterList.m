function varargout = FilterList(PredicateOrMask,varargin)
    
    Nfiltered = length(varargin);
    
    if (nargout ~= Nfiltered) && (nargout ~= 1)
        error('WWPP:Toolbox:InputOutputImbalance',...
            'The number of outputs should equal the number of vectors to be filtered.');
    end
    
    filtered = cellfun( @(v) v(PredicateOrMask,:),varargin,'UniformOutput',false);

    if (nargout == 1)
        varargout{1} = filtered;
    else
        varargout = filtered;
    end

end