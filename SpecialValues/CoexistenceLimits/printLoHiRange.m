function [] = printLoHiRange(range,name,performSort)
    
    if (nargin < 3)
        performSort = true;
    end
    
    if performSort
        range = sort(range);
    end
    
    if ischar(name)
        name = @(k) name ;
    elseif iscellstr(name)
        name = @(k) name{k};
    end
    
    fprintf('Smallest %s: %+23.16E\n',name(1),range(1));
    for k = 2:numel(range)-1
        fprintf('         %s: %+23.16E\n',name(k),range(k));
    end
    fprintf('Largest  %s: %+23.16E\n',name(numel(range)),range(end));
    

end