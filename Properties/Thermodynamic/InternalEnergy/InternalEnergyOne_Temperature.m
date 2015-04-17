function IntEnergy_T   = InternalEnergyOne_Temperature(rho,T)
    
    rhoc = CriticalDensity()        ;
    Tc   = CriticalTemperature()    ;
    
    delta = rho  / rhoc ;
    tau   = Tc  ./ T    ;
    
    IntEnergy_tau = InternalEnergyOneR_tau(delta,tau)   ;
    tau_T         = tau_Temperature(T)                  ;
    IntEnergy_T   = IntEnergy_tau .* tau_T              ;
    
end