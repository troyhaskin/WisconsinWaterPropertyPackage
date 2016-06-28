function [xSol,rSol,cSol] = NewtonUpdater(Update,Guess,Tolerance,MaxIter,~)
    
    % ======================================================================= %
    %                                   Set-Up                                %
    % ======================================================================= %
    xSol         = 0 * Guess             ;
    rSol         = xSol(:,1)             ;
    cSol         = xSol(:,1)             ;
    xk           = Guess                 ;
    N            = size(Guess,1)         ;
%     converged    = false(N,1)            ;
    Iupdate      = (1:N)'                ;
    SumErr       = xSol                  ;
    Iter         = 0                     ;
    alpha        = 0.5                   ; % Back-tracking relaxor
    
    
    % ======================================================================= %
    %     Pre-loop: check for already converged solution (just in case)       %
    % ======================================================================= %
    
    %   Account for a function handle that may return only one array
    try
        [dxNow,rBest] = Update(xk,Iupdate)    ; % First update and norm
    catch
        Update        = @(x,mask) functionHandleWorkAround(Update,x,mask)   ;
        [dxNow,rBest] = Update(xk,Iupdate)                                  ;
    end

    
    % Convergence check: Residual ONLY.
    converged       = abs(rBest) < Tolerance;
    notConverged    = not(converged)        ;
    notDone         = any(notConverged)     ;
    Ipush           = Iupdate(converged)	;
    xSol(Ipush,:)   = xk(converged,:)       ;
    rSol(Ipush)     = rBest(converged)      ;
    cSol(Ipush)     = true                  ;

    % Contract the unconverged values
    rBest       = rBest(notConverged)       ;
    Iupdate     = Iupdate(notConverged)     ;
    xk          = xk(notConverged,:)        ;
    dxNow       = dxNow(notConverged,:)     ;
    SumErr      = SumErr(notConverged,:)    ;
    
    
    
    % ======================================================================= %
    %                                Main Iterator                            %
    % ======================================================================= %
    while notDone
        
        % Update the system with full Newton step
        [dxNext,Rnew] = Update(xk - dxNow,Iupdate) ;
        
        % Back-track if needed by the following criteria
        NeedBackTrack = (rBest < Rnew)       | isnan(Rnew)              |...
                        any(isnan(dxNext),2) | any(imag(dxNext)~=0,2)   |...
                        any(imag(Rnew)~=0,2)                            ;
        
        % Back-tracker loop
        if any(NeedBackTrack)
            g = FilterList(NeedBackTrack,xk,dxNow,Iupdate,rBest);
            [dxNowBT,dxNextBT,RnewBT] = Backtrack(g{:},Update,alpha);
            dxNow (NeedBackTrack,:) = dxNowBT ;
            dxNext(NeedBackTrack,:) = dxNextBT;
            Rnew  (NeedBackTrack,:) = RnewBT  ;
        end
        
        % Set new values
        rBest = Rnew        ;
        xkp1  = xk - dxNow  ;
        
        %   Test for convergence
        belowTol  = abs(rBest)             < Tolerance          ;
        wontMove  = abs(sum(abs(dxNow),2)) < eps()              ;
        isNaN     = any(isnan([dxNow,rBest]),2)                 ;
        converged = belowTol | isNaN | isinf(rBest) | wontMove  ;
        
        %   Non-convergence stuff
        notConverged    = not(converged)                      ;
        BelowIterMax    = Iter < MaxIter                      ;
        Iter            = Iter + 1                            ;
        notDone         = any(notConverged) && BelowIterMax   ;
        
        % Push converged values into the solution vector
        Ipush                   = Iupdate(converged);
        xSol(Ipush,:)           = xkp1(converged,:) ;
        rSol(Ipush)             = rBest(converged)  ;
        cSol(Iupdate(belowTol)) = true              ;
        
        % Contract the unconverged values
        dxNow       = dxNext                    ;
        rBest       = rBest  (notConverged)     ;
        Iupdate     = Iupdate(notConverged)     ;
        xk          = xkp1   (notConverged,:)   ;
        dxNow       = dxNow  (notConverged,:)   ;
        SumErr      = SumErr (notConverged,:)   ;
    end
    
    xSol(Iupdate,:) = xk(:,:)   ;
    rSol(Iupdate)   = rBest     ;

end

function [dxNow,dxNext,Rnew] = Backtrack(xk,dxNow,iUpdate,Rbest,Update,alpha)
    
    [m,n]       = size(xk)  ;
    dxNext(m,n) = 0         ;
    Rnew(m,1)   = 0         ;
    relax       = true(m,1) ;
    
    while any(relax) && any(abs(dxNow(:)) > eps())
        dxNow(relax,:) = alpha * dxNow(relax,:);
        
        [dxNext(relax,:),Rnew(relax)] = Update(xk(relax,:) - dxNow(relax,:),iUpdate(relax))  ;
        
        relax(relax) =  (Rbest(relax) < Rnew(relax))    |...
                        isnan(Rnew(relax))              |...
                        any(isnan(dxNext(relax)),2)     |...
                        any(imag(dxNext(relax))~=0,2)   |...
                        any(imag(Rnew(relax))~=0,2)     ;
    end
end

function [dx,rNorm] = functionHandleWorkAround(f,x,mask)
    output = f(x,mask);
    dx = output(:,1:end-1);
    rNorm = output(:,end);
end