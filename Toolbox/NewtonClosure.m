function [xSol,Results] = NewtonClosure(Updater,Guess,Tolerance,MaxIter)
    
    % ======================================================================= %
    %                                   Set-Up                                %
    % ======================================================================= %
    
    % Initialize
    xSol    = 0 * Guess     ;
    xk      = Guess         ;
    N       = size(Guess,1) ;
    Iupdate = (1:N)'        ;
    SumErr  = xSol          ;
    Iter    = 0             ;
    alpha   = 0.1           ; % Back-tracking relaxor
    
    % Pull residual and residual derivative handles from the Updater
    
    
    
    
    % ======================================================================= %
    %     Pre-loop: check for already converged solution (just in case)       %
    % ======================================================================= %
    
    % Evaluate
    Rbest = Updater.GetResidual(xk);
    
    % Convergence check: Residual ONLY.
    Converged       = abs(Rbest) < Tolerance;
    NotConverged    = not(Converged)        ;
    NotDone         = any(NotConverged)     ;
    Ipush           = Iupdate(Converged)	;
    Rbest           = Rbest(NotConverged)   ;
    xSol(Ipush,:)   = xk(Converged,:)       ;
    
    % Contract the unconverged values
    Iupdate     = Iupdate(NotConverged)     ;
    xk          = xk(NotConverged,:)        ;
    SumErr      = SumErr(NotConverged,:)    ;
    
    
    
    % ======================================================================= %
    %                                Main Iterator                            %
    % ======================================================================= %
    while NotDone
        
        % Calculate the full Newton step
        Updater.SetFilter(Iupdate)          ;
        dRbest = Updater.GetDResidual(xk)   ;
        dx  = Rbest./dRbest                 ;
        
        % Get the new residual from the full step
        Rnew = Updater.GetResidual(xk - dx) ;
        
        % Back-track if needed by the following criteria
        NeedBackTrack = (abs(Rbest) < abs(Rnew)) | isnan(Rnew);
        
        % Back-tracker loop
        if any(NeedBackTrack)
            g = FilterList(NeedBackTrack,xk,dx,Iupdate,Rbest);
            [dx,Rnew] = AssignWithFilter(@() BackTracker(g{:},alpha,Updater),...
                NeedBackTrack,dx,Rnew);
            
            % Restore the closure's pre-back track filter
            Updater.SetFilter(Iupdate);
        end
        
        % Set new values
        Rbest = Rnew    ;
        xkp1  = xk - dx ;
        
        % Post-update loop-breaking checks
        Converged       = ConvergenceTest(dx,Rbest,Tolerance) ;
        NotConverged    = not(Converged)                      ;
        BelowIterMax    = Iter < MaxIter                      ;
        Iter            = Iter + 1                            ;
        NotDone         = any(NotConverged) && BelowIterMax   ;
        
        % Push converged values into the solution vector
        Ipush         = Iupdate(Converged)  ;
        xSol(Ipush,:) = xkp1(Converged,:)       ;
        
        % Contract the unconverged values
        Iupdate = Iupdate(NotConverged)     ;
        xk      = xkp1(NotConverged,:)      ;
        Rbest   = Rbest(NotConverged,:)     ;
        SumErr  = SumErr(NotConverged,:)    ;
    end
    
    Results = Updater.Finalize(xSol);
    
end

function [dx,Rnew] = BackTracker(xk,dx,Iupdate,Rbest,alpha,Updater)
    
    % Update the closure's filter
    Updater.SetFilter(Iupdate);
    
    % First relaxation
    dx            = alpha * dx                                  ;
    Rnew          = Updater.GetResidual(xk - dx)                ;
    NeedBackTrack = (abs(Rbest) < abs(Rnew))& (xk ~= (xk - dx)) ;
    
    % Back-tracking loop (via recursion)
    if any(NeedBackTrack)
        g = FilterList(NeedBackTrack,xk,dx,Iupdate,Rbest);
        [dx,Rnew] = AssignWithFilter(@() BackTracker(g{:},alpha,Updater),...
            NeedBackTrack,dx,Rnew);
    end
end

function Converged = ConvergenceTest(dx,Norm,Tolerance)
    IsZero      = abs(Norm)      < Tolerance        ;
    WontMove    = abs(sum(dx,2)) < Tolerance        ;
    IsNaN       = isnan(dx)      | isnan(Norm)      ;
    IsInf       = not(isfinite(Norm))               ;
    Converged	= IsZero | IsNaN | IsInf | WontMove ;
end
