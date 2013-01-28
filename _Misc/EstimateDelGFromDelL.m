function delG = EstimateDelGFromDelL(delL,tau0)
    
    Tolerance = 1E-13   ;
    IterMax   = 1E3     ;
    
    if (nargin < 2)
        tau0    = EstimateTauFromDelL(delL);
    end
    
    Updater = @(tau,Mask) GetTauFromDelL(tau,Mask,delL);
    tau     = NewtonUpdater(Updater,tau0,Tolerance,IterMax);
        
    delG = EstimateDelGFromTau(tau);
    
end

function [dx,Norm] = GetTauFromDelL(tau,Mask,delL)
    
    Residual  = EstimateDelLFromTau(tau) - delL(Mask);
    dResidual = EstimateDelLFromTau_tau(tau);
    
    dx   = Residual ./ dResidual;
    Norm = abs(Residual)        ;
    
end

