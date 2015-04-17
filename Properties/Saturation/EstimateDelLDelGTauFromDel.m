function [delL,delG,tau] = EstimateDelLDelGTauFromDel(del)
    
    % Allocate
    delL = del * 0;
    delG = del * 0;
    tau  = del * 0;
    
    % Phase logicals
    IsDelC   = (del == 1)   ;   % Critical condition
    IsDelL   = del >  1     ;   % liquid side in stable Newton iteration regime
    IsDelG   = del <  1     ;   % gas    side in stable Newton iteration regime
    
    % Push known values
    delL(IsDelL) = del(IsDelL);
    delG(IsDelG) = del(IsDelG);
    
    % Catch critical states and auto-assign (also avoids need for iteration)
    delL(IsDelC) = 1    ;
    delG(IsDelC) = 1    ;
    tau (IsDelC) = 1    ;
    
    % Iteration parameters
    Tolerance = 1E-10   ;   % Just above the FP precision limit of doubles
    IterMax   = DefaultMaximumIterationCount()   ;
    
    % Solve with liquid density
    if any(IsDelL)
        tauL = EstimateTauFromDelL(del(IsDelL)) ;   % Get initial guess
        
        Updater = @(t,Mask) GetTauFromDelL(t,Mask,del(IsDelL))  ;   % Update handle
        tauL = NewtonUpdater(Updater,tauL,Tolerance,IterMax)    ;   % Solve
        
        tau (IsDelL) = tauL                     ; % Assign tau value
        delG(IsDelL) = EstimateDelGFromTau(tauL); % Assign delG value
    end



    % Solve with gas density
    if any(IsDelG)
        tauG = EstimateTauFromDelG(del(IsDelG));    % Get initial guess
        
        Updater = @(t,Mask) GetTauFromDelG(t,Mask,del(IsDelG))  ;	% Update handle
        tauG = NewtonUpdater(Updater,tauG,Tolerance,IterMax)    ;   % Solve
        
        tau (IsDelG) = tauG                     ; % Assign tau value
        delL(IsDelG) = EstimateDelLFromTau(tauG); % Assign delL value
    end
    
end



function [dx,Norm] = GetTauFromDelL(tau,Mask,delL)
    
    Residual  = EstimateDelLFromTau(tau) - delL(Mask);
    dResidual = EstimateDelLFromTau_tau(tau);
    
    dx      = Residual ./ dResidual ;
    Norm = abs(Residual);
end

function [dx,Norm] = GetTauFromDelG(tau,Mask,delG)
    
    Residual  = EstimateDelGFromTau(tau) - delG(Mask);
    dResidual = EstimateDelGFromTau_tau(tau);
    
    dx   = Residual ./ dResidual;
    Norm = abs(Residual)        ;
    
end
