function dPdtau = Pressure_tau(rho,T,PhaseCheck)
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    %               [J/kg/K]                 [kg/m^3]              [K]  
    RrhoTc    = SpecificGasConstant() * CriticalDensity() * CriticalTemperature(); % [Pa]
            
    OnePhiHandle = @(delta,tau,Mask)(delta.*tau.*HelmholtzResidual_dt(delta,tau) - ...
                                     delta.*     HelmholtzResidual_d (delta,tau) -... 
                                     1).* delta ./ tau * RrhoTc;
                                 
    dPdtau = GenericTwoPhiProperty(rho,T,OnePhiHandle,'quality',PhaseCheck);
    
end