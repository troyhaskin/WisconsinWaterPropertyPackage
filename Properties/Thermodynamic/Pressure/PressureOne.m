function P = PressureOne(rho,T)
    
    rhoc = CriticalDensity();
    Tc   = CriticalTemperature();
    
    delta = rho  / rhoc;
    tau   = Tc  ./ T;
    
    P = PressureOneR(delta,tau);
    
end