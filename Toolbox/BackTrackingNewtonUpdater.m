function xSol = BackTrackingNewtonUpdater(Update,Objective,Guess,Tolerance,MaxIter,...
                                                                    Contraction,Relaxor)
    
    if (nargin < 7)
        Relaxor = 0.5;
    end
    
    if (nargin < 6)
        Contraction = true;
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
    
    switch(Contraction)
        case true
            SolveWithContraction();
            
        case false
            SolveWithoutContraction();
            
        otherwise
            error(['The optional input ''LoopContraction'' must be ',...
                'a logical scalar (true/false).']);
    end
    
    
    function [] = SolveWithContraction()
        
        [dx,RNorm] = Update(xk,Iupdate); % Initial search direction and RNorm

        while NotDone
            
            %       Backtracking
            % =======================
            
            % Form the candidate solution
            [xkp1,SumErr] = KahanSum(xk,-dx,SumErr) ;
            alpha          = 1                      ;
            StepSizeTooBig = true                   ;
            
            % Back-tracking/Newton-relaxing loop
            while StepSizeTooBig
                RNormkp1        = Objective(xkp1,Iupdate)               ;
                ErrorNotReduced = (RNormkp1 > RNorm)                    ;
                StepSizeTooBig  = ErrorNotReduced || isnan(RNormkp1)    ;

                if StepSizeTooBig
                    alpha         = Relaxor * alpha                 ;
                    [xkp1,SumErr] = KahanSum(xk,-alpha*dx,SumErr)   ;
                    disp(['Backtracking:: alpha: ',num2str(alpha,'%+15.8E'),'   ',...
                                          'Iteration: ',num2str(Iter,'%+04G')]);
                end
            end
            
            % Update search direction and RNorm
            [dx,RNorm] = Update(xkp1,Iupdate);



            %       Convergence Checking and Vector Contraction
            % =======================================================
            
            % Post-update loop-breaking checks
            Converged       = ConvergenceTest(alpha*dx./xkp1,RNorm,Tolerance) ;
            NotConverged    = not(Converged)                            ;
            BelowIterMax    = Iter < MaxIter                            ;
            Iter            = Iter + 1                                  ;
            NotDone         = any(NotConverged) && BelowIterMax         ;

            % Push converged values into the solution vector
            Ipush         = Iupdate(Converged)  ;
            xSol(Ipush,:) = xkp1(Converged,:)   ;

            % Contract the unconverged values
            Iupdate     = Iupdate(NotConverged  )   ;
            xk          = xkp1   (NotConverged,:)   ;
            SumErr      = SumErr (NotConverged,:)   ;
        end
    end
    
    
    function [] = SolveWithoutContraction()
        
        [dx,RNorm] = Update(xk); % Initial search direction and RNorm
        
        while NotDone
            
            %       Backtracking
            % =======================
            
            % Form the candidate solution
            [xkp1,SumErr] = KahanSum(xk,-dx,SumErr) ;
            alpha          = 1                      ;
            StepSizeTooBig = true                   ;
            
            % Back-tracking/Newton-relaxing loop
            while StepSizeTooBig
                RNormkp1        = Objective(xkp1)                                   ;
                ErrorNotReduced = (RNormkp1 > RNorm)                                ;
                StepSizeTooBig  = ErrorNotReduced || isnan(NormResidualCandidate)   ;

                if StepSizeTooBig
                    alpha         = Relaxor * alpha                 ;
                    [xkp1,SumErr] = KahanSum(xk,-alpha*dx,SumErr)   ;
                end
            end
            
            % Update search direction and RNorm
            [dx,RNorm] = Update(xkp1);



            %                  Convergence Checking
            % =======================================================

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
