function cv = HeatCapacityIsochoricOne(rho,T)
    
    rhoc = CriticalDensity();
    Tc   = CriticalTemperature();
    
    delta = rho  / rhoc;
    tau   = Tc  ./ T;
    
    cv = HeatCapacityIsochoricOneR(delta,tau);
    
end