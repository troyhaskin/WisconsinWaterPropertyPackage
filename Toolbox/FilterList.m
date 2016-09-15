function varargout = FilterList(PredicateOrMask,varargin)
    
    Nfiltered = length(varargin);
    
    if (nargout ~= Nfiltered) && (nargout ~= 1)
        error('WWPP:Toolbox:InputOutputImbalance',...
            'The number of outputs should equal the number of vectors to be filtered.');
    end
    
    %   Apply the filter with scalar expansion through sift.
    filtered = cellfun( @(v) sift(v,PredicateOrMask),varargin,'UniformOutput',false);

    if (nargout == 1)
        varargout{1} = filtered;
    else
        varargout = filtered;
    end

end

function v = sift(v,mask)
    if isscalar(v)
        v = v*size(numel(mask),1);
    else
        v = v(mask,:);
    end
end