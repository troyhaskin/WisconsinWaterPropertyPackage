function [delL,delG,tau] = EstimateDelLDelGTauFromDel(del)
    
    % Allocate
    delL = del * 0;
    delG = del * 0;
    tau  = del * 0;
    
    % Phase logicals
    IsDelC = del == 1; % critical state
    IsDelL = del >  1; % liquid side
    IsDelG = del <  1; % gas    side
    
    % Push known values
    delL(IsDelL) = del(IsDelL);
    delG(IsDelG) = del(IsDelG);
    
    % Catch critical states and auto-assign (also avoids need for iteration)
	delL(IsDelC) = 1    ;
    delG(IsDelC) = 1    ;
    tau (IsDelC) = 1    ;
    
    % Iteration parameters
    Tolerance = 1E-13   ;
    IterMax   = 1E3     ;
    
    % Solve with liquid density
    if any(IsDelL)
        tauL0       = EstimateTauFromDelL(del(IsDelL));                 % Get initial guess
        Updater     = @(tau,Mask) GetTauFromDelL(tau,Mask,del(IsDelL)); % Update handle
        tauL        = NewtonUpdater(Updater,tauL0,Tolerance,IterMax);   % Solve
        delGsol     = EstimateDelGFromTau(tauL);                        % Back substitute
        
        tau (IsDelL) = tauL    ; % Assign tau value
        delG(IsDelL) = delGsol ; % Assign delG value
    end
    
    % Solve with gas density
    if any(IsDelG)
        tauG0       = EstimateTauFromDelG(del(IsDelG));                 % Get initial guess
        Updater     = @(tau,Mask) GetTauFromDelG(tau,Mask,del(IsDelG)); % Update handle
        tauG        = NewtonUpdater(Updater,tauG0,Tolerance,IterMax);   % Solve
        delLsol     = EstimateDelLFromTau(tauG);                        % Back substitute
        
        tau (IsDelG) = tauG    ; % Assign tau value
        delL(IsDelG) = delLsol ; % Assign delL value
    end

end

function [dx,Norm] = GetTauFromDelL(tau,Mask,delL)
    
    Residual  = EstimateDelLFromTau(tau) - delL(Mask);
    dResidual = EstimateDelLFromTau_tau(tau);

    dx      = Residual ./ dResidual ;
    Norm = abs(dx);
end

function [dx,Norm] = GetTauFromDelG(tau,Mask,delG)
    
    Residual  = EstimateDelGFromTau(tau) - delG(Mask);
    dResidual = EstimateDelGFromTau_tau(tau);
    
    dx   = Residual ./ dResidual;
    Norm = abs(Residual)        ;
    
end


