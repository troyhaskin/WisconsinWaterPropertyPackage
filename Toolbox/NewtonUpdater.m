function xSol = NewtonUpdater(Update,Guess,Tolerance,MaxIter,~)
    
    if (nargin < 5)
%         LoopContraction = true;
    end
    
    xSol            = 0 * Guess             ;
    xk              = Guess                 ;
    N               = size(Guess,1)         ;
    Converged       = false(N,1)            ;
    NotConverged    = not(Converged)        ;
    Iupdate         = (1:N)'                ;
    Iupdate         = Iupdate(NotConverged) ;
    SumErr          = xSol                  ;
    Iter            = 0                     ;
    NotDone         = true                  ;
    
    
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

function Converged = ConvergenceTest(dx,Norm,Tolerance)
    IsZero      = abs(Norm)   < Tolerance           ;
    WontMove    = any(abs(sum(dx,2)) < Tolerance)   ;
    IsNaN       = any(isnan(dx),2) | isnan(Norm)    ;
    IsInf       = not(isfinite(Norm))               ;
    Converged	= IsZero | IsNaN | IsInf | WontMove ;
end
