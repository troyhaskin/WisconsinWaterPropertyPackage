function [Uniques,iUniques,iVector] = UniqueEnough(Vector,AbsoluteTolerance)

    if (nargin == 1) || isempty(AbsoluteTolerance)
        AbsoluteTolerance = eps();
    end

    %   Set-Up and allocation
    Uniques  = Vector           ;
    N        = length(Vector)   ;
    iUniques = (1:N)'           ;
    iVector  = iUniques         ;

    %   Loop Set-Up
    k       = 1     ; % Unique counter
    NotDone = true  ; % Loop breaker

    while NotDone
        
        %   Grab the next non-unique enough values
        TestValue = Uniques(k);
        
        %   Form the logical comparison mask against remaining uniques
        IsUniqueEnough = abs(TestValue - Uniques) <= AbsoluteTolerance;


        %   Tests for infinites (Inf is considered unique of -Inf) and not-a-numbers
        if none(IsUniqueEnough)
            if isnan(TestValue)
                IsUniqueEnough = isnan(Uniques);
            elseif isinf(TestValue)
                
                if (TestValue > 0)
                    IsUniqueEnough = isinf(Uniques) & (Uniques > 0);
                else
                    IsUniqueEnough = isinf(Uniques) & (Uniques < 0);
                end
            end
        end

        
        %   Assign the current unique count to all appropriate places in the
        %   index vector which regenerates Vector when used as an index in Uniques
        iVector(iUniques(IsUniqueEnough)) = k;


        %   This ensures that the current unique index is not contracted out.
        IsUniqueEnough(k) = false;
  
      
        %   Contract the index vector which generates Uniques when used as an index for Vector
        iUniques = iUniques(not(IsUniqueEnough));
 
       
        %   Reassign the Uniques vector with current unique index list
        %   This is essentially a contraction if multiple uniques were found this loop.
        Uniques = Vector(iUniques);


        %   If the current index of Vector being checked is the last index, the 
        %   loop is ended. Otherwise, the counter is incremented and we iterate.
        if (iUniques(k) < N) && (k < length(Uniques))
            k = k + 1;
        else
            NotDone = false;
        end

    end

end