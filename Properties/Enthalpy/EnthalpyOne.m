function h = EnthalpyOne(rho,T)
    
    rhoc = CriticalDensity();
    Tc   = CriticalTemperature();
    
    delta = rho  / rhoc;
    tau   = Tc  ./ T   ; 
    
    h = EnthalpyOneR(delta,tau);
    
end