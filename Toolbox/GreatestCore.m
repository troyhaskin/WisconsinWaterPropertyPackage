function Greatest = GreatestCore(IsGreater,varargin)
    
    Greatest = varargin{1}      ;
    Nin      = length(varargin) ;
    
    for  k = 2:Nin
        Candidate = varargin{k};
        if IsGreater(Candidate,Greatest)
            Greatest = Candidate;
        end
    end
    
end