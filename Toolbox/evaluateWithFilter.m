function output = evaluateWithFilter(functionCall,filter,varargin)
    
    for k = 1:length(varargin)
        varargin{k} = varargin{k}(filter,:);
    end
    
    output = functionCall(varargin{:});
    
end
    