function w = SoundSpeedOne(rho,T)
    
    rhoc = CriticalDensity();
    Tc   = CriticalTemperature();
    
    delta = rho  / rhoc;
    tau   = Tc  ./ T;
    
    w = SoundSpeedOneR(delta,tau);
    
end