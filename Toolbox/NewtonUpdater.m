function xSol = NewtonUpdater(Update,Guess,Tolerance,MaxIter,~)
    
    % ======================================================================= %
    %                                   Set-Up                                %
    % ======================================================================= %
    xSol            = 0 * Guess             ;
    xk              = Guess                 ;
    N               = size(Guess,1)         ;
    Converged       = false(N,1)            ;
    NotConverged    = not(Converged)        ;
    Iupdate         = (1:N)'                ;
    Iupdate         = Iupdate(NotConverged) ;
    SumErr          = xSol                  ;
    Iter            = 0                     ;
    alpha           = 0.5                   ; % Back-tracking relaxor
    
    
    
    % ======================================================================= %
    %     Pre-loop: check for already converged solution (just in case)       %
    % ======================================================================= %
    
    % Evaluate
    [dxNow,Rbest] = Update(xk,Iupdate)    ; % First update and norm
    
    % Convergence check: Residual ONLY.
    Converged       = abs(Rbest) < Tolerance;
    NotConverged    = not(Converged)        ;
    NotDone         = any(NotConverged)     ;
    Ipush           = Iupdate(Converged)	;
    xSol(Ipush,:)   = xk(Converged,:)       ;
    
    % Contract the unconverged values
    Iupdate     = Iupdate(NotConverged)     ;
    xk          = xk(NotConverged,:)        ;
    SumErr      = SumErr(NotConverged,:)    ;
    
    
    
    % ======================================================================= %
    %                                Main Iterator                            %
    % ======================================================================= %
    while NotDone
        
        % Update the system with full Newton step
        [dxNext,Rnew] = Update(xk - dxNow,Iupdate) ;
        
        % Back-track if needed by the following criteria
        NeedBackTrack = (Rbest < Rnew) | isnan(Rnew) | any(isnan(dxNext),2) | ...
                              any(imag(dxNext)~=0,2) | any(imag(Rnew)~=0,2);
        
        % Back-tracker loop
        if any(NeedBackTrack)
            g = FilterList(NeedBackTrack,xk,dxNow,Iupdate,Rbest);
            [dxNow,dxNext,Rnew] = AssignWithFilter(@() BackTracker(g{:},alpha,Update),...
                NeedBackTrack,dxNow,dxNext,Rnew);
        end
        
        % Set new values
        Rbest = Rnew        ;
        xkp1  = xk - dxNow  ;
        dxNow = dxNext      ;
        
        % Post-update loop-breaking checks
        Converged       = ConvergenceTest(dxNow,Rbest,Tolerance) ;
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

function [dxNow,dxNext,Rnew] = BackTracker(xk,dxNow,Iupdate,Rbest,alpha,Update)
    
    % First relaxation
    dxNow         = alpha * dxNow               ;
    [dxNext,Rnew]  = Update(xk - dxNow,Iupdate) ;
    NeedBackTrack = (Rbest < Rnew) | isnan(Rnew) | any(isnan(dxNext),2) | ...
                    any(imag(dxNext)~=0,2) | any(imag(Rnew)~=0,2);
    
    % Back-tracking loop (via recursion)
    if any(NeedBackTrack)
        g = FilterList(NeedBackTrack,xk,dxNow,Iupdate,Rbest);
        [dxNow,dxNext,Rnew] = AssignWithFilter(@() BackTracker(g{:},alpha,Update),...
            NeedBackTrack,dxNow,dxNext,Rnew);
    end
end

function Converged = ConvergenceTest(dx,Norm,Tolerance)
    IsZeroAbs   = abs(Norm)    < Tolerance          ;
    WontMove    = any(abs(sum(dx,2)) < Tolerance)   ;
    IsNaN       = any(isnan(dx),2) | isnan(Norm)    ;
    IsInf       = not(isfinite(Norm))               ;
    Converged	= IsZeroAbs | IsNaN | IsInf | WontMove ;
end
