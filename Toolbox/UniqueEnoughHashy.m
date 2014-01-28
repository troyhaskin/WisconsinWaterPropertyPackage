function [Uniques,iUniques,iVector] = UniqueEnoughHashy(Vector,AbsoluteTolerance)

    if (nargin == 1) || isempty(AbsoluteTolerance)
%         AbsoluteTolerance = eps();
    end

    %   Set-Up and allocation
    Uniques  = Vector           ;
    N        = length(Vector)   ;
    iUniques = (1:N)'           ;
    iVector  = iUniques         ;
    
    UniqueCount = 0;
    StringNumbers = char(regexprep(cellstr(num2str(Vector(:),'%22.16E')),'(\.|E[\-\+])',''));
    Hash = struct();
    
    for k = 1:N
        HashKey = ['x',StringNumbers(k,:)];
        
        if isfield(Hash,HashKey)
            iVector(k) = Hash.(HashKey);
        else
            UniqueCount           = UniqueCount + 1;
            iVector(k)            = UniqueCount;
            Uniques(UniqueCount)  = Vector(k);
            iUniques(UniqueCount) = k;
            Hash.(HashKey)        = UniqueCount;
        end
    end
    
    Uniques  =  Uniques(1:UniqueCount);
    iUniques = iUniques(1:UniqueCount);

end