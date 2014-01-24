function varargout = AssignWithFilter(FunctionCall,Filter,varargin)
    
    % Allocation
    N         = length(varargin);
    varargout = varargin;
    
    % Function-return pattern
    switch(N)
        case 0
            varargout = cell(1,0);
        case 1
            x1 = FunctionCall();
            varargout{1}(Filter,:) = x1;
            
        case 2
            [x1,x2] = FunctionCall();
            varargout{1}(Filter,:) = x1;
            varargout{2}(Filter,:) = x2;
            
        case 3
            [x1,x2,x3] = FunctionCall();
            varargout{1}(Filter,:) = x1;
            varargout{2}(Filter,:) = x2;
            varargout{3}(Filter,:) = x3;
            
        case 4
            [x1,x2,x3,x4] = FunctionCall();
            varargout{1}(Filter,:) = x1;
            varargout{2}(Filter,:) = x2;
            varargout{3}(Filter,:) = x3;
            varargout{4}(Filter,:) = x4;
            
        case 5
            [x1,x2,x3,x4,x5] = FunctionCall();
            varargout{1}(Filter,:) = x1;
            varargout{2}(Filter,:) = x2;
            varargout{3}(Filter,:) = x3;
            varargout{4}(Filter,:) = x4;
            varargout{5}(Filter,:) = x5;
        case 6
            [x1,x2,x3,x4,x5,x6] = FunctionCall();
            varargout{1}(Filter,:) = x1;
            varargout{2}(Filter,:) = x2;
            varargout{3}(Filter,:) = x3;
            varargout{4}(Filter,:) = x4;
            varargout{5}(Filter,:) = x5;
            varargout{6}(Filter,:) = x6;
    end
    
    
end