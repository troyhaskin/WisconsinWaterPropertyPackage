function [xSol,Results] = NewtonClosure(ClosureConstructor,Guess,Tolerance,MaxIter)
    
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
    
    % Create closures
    UpdateMain = ClosureConstructor();  % Used in the outer loop for Newton updating
    UpdateBack = ClosureConstructor();  % Used in the inner loop for back-tracking



    % ======================================================================= %
    %     Pre-loop: check for already converged solution (just in case)       %
    % ======================================================================= %
    
    % Loop vectors allocation
    Converged       = false(N,1)            ;
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
        
        % Calculate the full Newton step
        UpdateMain.SetFilter(Iupdate)       ;
        Rbest  = UpdateMain.GetResidual(xk) ;
        dRbest = UpdateMain.GetDResidual(xk);
        dx     = Rbest./dRbest              ;
        
        % Get the new residual from the full step
        Rnew = UpdateMain.GetResidual(xk - dx) ;
        
        % Back-track if needed by the following criteria
        NeedBackTrack = (abs(Rbest) < abs(Rnew)) | isnan(Rnew);
        
        % Back-tracker loop
        if any(NeedBackTrack)
            g  = FilterList(NeedBackTrack,xk,dx,Iupdate,Rbest);
            dx = AssignWithFilter(@() BackTracker(g{:},alpha,UpdateBack),...
                NeedBackTrack,dx);
        end
        
        % Set new values
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
        SumErr  = SumErr(NotConverged,:)    ;
    end
    
    Results = UpdateMain.Finalize(xSol);
    
end

function dx = BackTracker(xk,dx,Iupdate,Rbest,alpha,Update)
    
    % Update the closure's filter
    Update.SetFilter(Iupdate);
    
    % First relaxation
    dx            = alpha * dx                                  ;
    Rnew          = Update.GetResidual(xk - dx)                ;
    NeedBackTrack = (abs(Rbest) < abs(Rnew))& (xk ~= (xk - dx)) ;
    
    % Back-tracking loop (via recursion)
    if any(NeedBackTrack)
        g  = FilterList(NeedBackTrack,xk,dx,Iupdate,Rbest);
        dx = AssignWithFilter(@() BackTracker(g{:},alpha,Update),...
            NeedBackTrack,dx);
    end
end

function Converged = ConvergenceTest(dx,Norm,Tolerance)
    IsZero      = abs(Norm)      < Tolerance        ;
    WontMove    = abs(sum(dx,2)) < Tolerance        ;
    IsNaN       = isnan(dx)      | isnan(Norm)      ;
    IsInf       = not(isfinite(Norm))               ;
    Converged	= IsZero | IsNaN | IsInf | WontMove ;
end
