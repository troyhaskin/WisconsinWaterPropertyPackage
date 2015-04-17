function T_tau = Temperature_tau(tau)
    
    Tc    = CriticalTemperature()   ;
    T_tau = -Tc ./ tau.^2           ;
    
end