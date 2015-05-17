function varargout = AssignWithFilter(functionCall,predicate,varargin)
    
    % Allocation
    n               = length(varargin)          ;
    functionOutputs = cell(1,length(varargin))  ;
    varargout       = varargin                  ;
    
    %   Call function assuming an output arity of n
    if (nargin(functionCall) ~= 0)
        [functionOutputs{:}] = functionCall(predicate);
    else
        [functionOutputs{:}] = functionCall();
    end
    
    %   Update varargout with function results according to the predicate
    for k = 1:n
        varargout{k}(predicate,:) = functionOutputs{k};
    end
    
end