function [delL,delG,tau] = EstimateDelLDelGTauFromDel(del,isDelL,isDelG)
    
    % Allocate
    delL = del * 0;
    delG = del * 0;
    tau  = del * 0;
    
    % Phase logicals
    isDelC   = (del == 1)   ;   % Critical condition
    if (nargin < 2)
        isDelL   = del > 1  ;   % liquid side in stable Newton iteration regime
        isDelG   = del < 1  ;   % gas    side in stable Newton iteration regime
    end
    
    % Push known values
    delL(isDelL) = del(isDelL);
    delG(isDelG) = del(isDelG);
    
    % Catch critical states and auto-assign (also avoids need for iteration)
    delL(isDelC) = 1    ;
    delG(isDelC) = 1    ;
    tau (isDelC) = 1    ;
    
    % Iteration parameters
    Tolerance = 1E-14                           ;   % Just above the FP precision limit of doubles
    IterMax   = DefaultMaximumIterationCount()  ;
    
    % Solve with liquid density
    if any(isDelL)
        tauL = EstimateTauFromDelL(del(isDelL)) ;   % Get initial guess
        
        Updater = @(tau,Mask) GetTauFromDelL(tau,Mask,del(isDelL)) ;   % Update handle
        tauL = NewtonUpdater(Updater,tauL,Tolerance,IterMax);   % Solve
        
        tau (isDelL) = tauL                     ; % Assign tau value
        delG(isDelL) = EstimateDelGFromTau(tauL); % Assign delG value
    end



    % Solve with gas density
    if any(isDelG)
        tauG = EstimateTauFromDelG(del(isDelG));    % Get initial guess
        
        Updater = @(t,Mask) GetTauFromDelG(t,Mask,del(isDelG))  ;   % Update handle
        tauG = NewtonUpdater(Updater,tauG,Tolerance,IterMax)    ;   % Solve
        
        tau (isDelG) = tauG                     ; % Assign tau value
        delL(isDelG) = EstimateDelLFromTau(tauG); % Assign delL value
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
