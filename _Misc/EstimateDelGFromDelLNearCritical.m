function [delL,delG] = EstimateDelGFromDelLNearCritical(del)
    
    % The NewtonUpdater assumes column vectors
    [del,SizeDel] = Columnify(del);
    
    % Allocate
    delL = del * 0;
    delG = del * 0;
    
    % Phase logicals
    IsDelL = del >= 1;
    IsDelG = del <= 1;
    IsDelC = del == 1;
    
    % Push known values
    delL(IsDelL) = del(IsDelL);
    delG(IsDelG) = del(IsDelG);
    
    % Initial tau guess that is used as an intermediary
    tauL0 = EstimateTauFromDelL(del(IsDelL));
    tauG0 = EstimateTauFromDelG(del(IsDelG));
    
    % Catch critical states and auto-assign for instant convergence
    tauL0(IsDelC) = 1	;
    tauG0(IsDelC) = 1	;
   
    % Iteration parameters
    Tolerance = 1E-13   ;
    IterMax   = 1E3     ;
    UpdaterL   = @(tau,Mask) GetTauFromDelL(tau,Mask,del(IsDelL));
    UpdaterG   = @(tau,Mask) GetTauFromDelG(tau,Mask,del(IsDelG));
    
    % Solve
    tauL = NewtonUpdater(UpdaterL,tauL0,Tolerance,IterMax);
    tauG = NewtonUpdater(UpdaterG,tauG0,Tolerance,IterMax);
        
    % Back substitute
    delLsol = EstimateDelLFromTau(tauL);
    delGsol = EstimateDelGFromTau(tauG);
    
    % Push solved values into output vectors
    delL(not(IsDelL)) = delLsol;
    delG(not(IsDelG)) = delGsol;
    
    % Restore input shape
    delL = RestoreShape(delL,SizeDel);
    delG = RestoreShape(delG,SizeDel);
end

function [dx,Norm] = GetTauFromDelL(tau,Mask,delL)
    
    Residual  = EstimateDelLFromTau(tau) - delL(Mask);
    dResidual = EstimateDelLFromTau_tau(tau);
    
    dx   = Residual ./ dResidual;
    Norm = abs(Residual)        ;
    
end

function [dx,Norm] = GetTauFromDelG(tau,Mask,delG)
    
    Residual  = EstimateDelGFromTau(tau) - delG(Mask);
    dResidual = EstimateDelGFromTau_tau(tau);
    
    dx   = Residual ./ dResidual;
    Norm = abs(Residual)        ;
    
end


