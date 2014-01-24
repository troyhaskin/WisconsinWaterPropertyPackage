function varargout = ApplyMask(Mask,varargin)
    N      = length(varargin);
    
    switch (isa(varargin{1},'function_handle'))
        case true
            Masked    = cellfun( @(v) v(Mask),varargin(2:end),'UniformOutput',false)    ;
            Callback  = varargin{1}                                                     ;
            varargout = Callback(Masked{:})                                             ;


        case false
            Masked = cellfun( @(v) v(Mask),varargin,'UniformOutput',false);
            
            if (nargout == 1)  || (nargout == 0)
                varargout{1} = Masked;
            end
            if (nargout == N)
                varargout = cell(1,N);
                
                for k = 1:N
                    varargout{k} = Masked{k};
                end
            end

    end
end

