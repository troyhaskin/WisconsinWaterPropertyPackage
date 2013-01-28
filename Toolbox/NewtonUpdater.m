function xSol = NewtonUpdater(Update,Guess,Tolerance,MaxIter,LoopContraction)
    
    if (nargin < 5)
        LoopContraction = true;
    end
    
    xSol            = 0 * Guess             ;
    xk              = Guess                 ;
    N               = size(Guess,1)         ;
    Converged       = false(1,N)            ;
    NotConverged    = not(Converged)        ;
    Iupdate         = 1:N                   ;
    Iupdate         = Iupdate(NotConverged) ;
    SumErr          = xSol                  ;
    Iter            = 0                     ;
    NotDone         = true                  ;
    
    switch(LoopContraction)
        case true
            SolveWithContraction();
            
        case false
            SolveWithoutContraction();
            
        otherwise
            error(['The optional input ''LoopContraction'' must be ',...
                'a logical scalar (true/false).']);
    end
    
    
    function [] = SolveWithContraction()
        
        while NotDone
            
            % Update the system
            [dx,RNorm] = Update(xk,Iupdate) ;
            [xkp1,SumErr] = KahanSum(xk,-dx,SumErr);
            
            % Post-update loop-breaking checks
            Converged       = ConvergenceTest(dx,RNorm,Tolerance) ;
            NotConverged    = not(Converged)                      ;
            BelowIterMax    = Iter < MaxIter                      ;
            Iter            = Iter + 1                            ;
            NotDone         = any(NotConverged) && BelowIterMax   ;
            
            % Push converged values into the solution vector
            Ipush         = Iupdate(Converged)  ;
            xSol(Ipush,:) = xkp1(Converged,:)       ;
            
            % Contract the unconverged values
            Iupdate     = Iupdate(NotConverged)     ;
            xk          = xkp1(NotConverged,:)      ;
            SumErr      = SumErr(NotConverged,:)    ;
        end
    end
    
    
    function [] = SolveWithoutContraction()
        
        while NotDone
            
            % Update the system
            [dx,RNorm] = Update(xk) ;
            [xkp1,SumErr] = KahanSum(xk,-dx,SumErr);
            
            % Post-update loop-breaking checks
            Converged       = ConvergenceTest(dx./xk,RNorm,Tolerance)   ;
            NotConverged    = not(Converged)                            ;
            BelowIterMax    = Iter < MaxIter                            ;
            Iter            = Iter + 1                                  ;
            NotDone         = any(NotConverged) && BelowIterMax         ;
            
            % Push converged values into the solution vector
            Ipush       = Iupdate(Converged)    ;
            xSol(Ipush) = xkp1(Ipush)           ;
        end
    end
end

function Converged = ConvergenceTest(dx,Norm,Tolerance)
    IsZero      = abs(Norm)   < Tolerance           ;
    WontMove    = any(abs(dx) < Tolerance,2)        ;
    IsNaN       = any(isnan(dx),2) | isnan(Norm)    ;
    IsInf       = not(isfinite(Norm))               ;
    Converged	= IsZero | IsNaN | IsInf | WontMove ;
end
