function P_delta = PressureSingle_delta(rho,T)
    
    rhoc = CriticalDensity()    ; %[kg/m^3]
    Tc   = CriticalTemperature(); %[K]
    
    delta = rho  / rhoc ;
    tau   = Tc  ./ T    ;
    
    P_delta = PressureSingleR_delta(delta,tau);
end