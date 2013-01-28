function dPdi = Pressure_InternalEnergy(rho,T,PhaseCheck)
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    Tc        = CriticalTemperature()   ; %[K]
    R         = SpecificGasConstant()   ; %[J/kg/K]        
            
    OnePhiHandle = @(delta,tau,Mask)(HelmholtzIdealGas_t(delta,tau)  +...
                                     HelmholtzResidual_t(delta,tau)) .*... 
                                     R .* Tc ;
                                 
    dPdi = GenericTwoPhiProperty(rho,T,OnePhiHandle,'quality',PhaseCheck);
    
end