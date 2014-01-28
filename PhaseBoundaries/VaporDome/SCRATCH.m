function Closure = SCRATCH()

    a = 1;
    
    Closure.set = @(x) set(x);
    Closure.get = @()  get();
    
    function [] = set(x)
        a = x;
    end
    
    function x = get()
        x = a;
    end

end
