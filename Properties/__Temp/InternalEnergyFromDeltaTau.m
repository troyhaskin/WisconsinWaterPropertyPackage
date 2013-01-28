function i = InternalEnergyFromDeltaTau(delta,tau,PhaseCheck)
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    Tc        = CriticalTemperature()   ; %[K]
    R         = SpecificGasConstant()   ; %[J/kg/K]        
            
    OnePhiHandle = @(delta,tau,Mask) SinglePhase_tau(delta,tau,R*Tc);                                
    i = GenericTwoPhiPropertyFromDeltaTau(delta,tau,OnePhiHandle,'quality',PhaseCheck);
    
end

function i = SinglePhase_tau(delta,tau,RTc)
    i = (HelmholtzIdealGas_t(delta,tau)    +...
         HelmholtzResidual_t(delta,tau))   .* RTc;
end

