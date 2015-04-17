function lam = ThermalConductivityOne(rho,T)
    
    rhoc = CriticalDensity();
    Tc   = CriticalTemperature();
    
    delta = rho  / rhoc;
    tau   = Tc  ./ T;
    
    lam = ThermalConductivityOneR(delta,tau);
end