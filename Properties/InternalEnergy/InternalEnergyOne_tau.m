function IntEnergy_t   = InternalEnergyOne_tau(rho,T)
    
    rhoc = CriticalDensity();
    Tc   = CriticalTemperature();
    
    delta = rho  / rhoc;
    tau   = Tc  ./ T;
    
    IntEnergy_t = InternalEnergyOneR_tau(delta,tau);
    
end