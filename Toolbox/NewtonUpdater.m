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
    Rbest       = Rbest(NotConverged)       ;
    Iupdate     = Iupdate(NotConverged)     ;
    xk          = xk(NotConverged,:)        ;
    dxNow       = dxNow(NotConverged,:)     ;
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
            [dxNowBT,dxNextBT,RnewBT] = Backtrack(g{:},Update,alpha);
            dxNow (NeedBackTrack) = dxNowBT ;
            dxNext(NeedBackTrack) = dxNextBT;
            Rnew  (NeedBackTrack) = RnewBT  ;
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
        Rbest       = Rbest  (NotConverged)     ;
        Iupdate     = Iupdate(NotConverged)     ;
        xk          = xkp1   (NotConverged,:)   ;
        dxNow       = dxNow  (NotConverged,:)   ;
        SumErr      = SumErr (NotConverged,:)   ;
    end
    
    Ipush         = Iupdate(NotConverged)  ;
    xSol(Ipush,:) = xk(NotConverged,:)     ;

end

function [dxNowBT,dxNextBT,RnewBT] = Backtrack(xk,dxNow,iUpdate,Rbest,Update,alpha)

    dxNowBT  = xk;
    dxNextBT = xk;
    RnewBT   = xk;
    iDone    = 1:length(xk);
    
    while not(isempty(iUpdate))
        dxNow         = alpha * dxNow               ;

        [dxNext,Rnew] = Update(xk - dxNow,iUpdate)  ;
        NeedBackTrack = (Rbest < Rnew) | isnan(Rnew) | any(isnan(dxNext),2) | ...
                        any(imag(dxNext)~=0,2) | any(imag(Rnew)~=0,2);
    
        iPush   = iDone(not(NeedBackTrack)) ;
        iUpdate = iUpdate(NeedBackTrack)    ;
        
        dxNowBT(iPush)  = dxNow(iPush);
        dxNextBT(iPush) = dxNext(iPush);
        RnewBT(iPush)   = Rnew(iPush);
    end
end

function Converged = ConvergenceTest(dx,Norm,Tolerance)
    IsZeroAbs   = abs(Norm)     < Tolerance             ;
    WontMove    = abs(sum(dx,2)) < Tolerance            ;
    IsNaN       = any(isnan(dx),2) | isnan(Norm)        ;
    IsInf       = not(isfinite(Norm))                   ;
    Converged	= IsZeroAbs | IsNaN | IsInf | WontMove  ;
end
