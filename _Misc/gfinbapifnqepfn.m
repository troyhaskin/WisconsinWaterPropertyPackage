function delG = EstimateDelGFromDelLNearCritical(delL)
    
    Tolerance = 1E-12   ;
    IterMax   = 1E3     ;
    
    
    tau0    = EstimateTauFromDelL(delL);
    Updater = @(tau,Mask) GetTauFromDelL(tau,Mask,delL);
    tau     = NewtonUpdater(Updater,tau0,Tolerance,IterMax);
        
    delG = EstimateDelGFromTau(tau);
    
end

function [dx,Norm] = GetTauFromDelL(tau,Mask,delL)
    
    Residual  = EstimateDeLLFromTau(tau) - delL(Mask);
    dResidual = EstimateDelLFromTau_tau(tau);
    
    dx   = Residual ./ dResidual;
    Norm = abs(Residual)        ;
    
end

