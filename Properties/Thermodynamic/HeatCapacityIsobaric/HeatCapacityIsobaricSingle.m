function cp = HeatCapacityIsobaricSingle(rho,T)
    % Single phase property with natural state specification
    
    Tc    = CriticalTemperature()   ; %[K]
    rhoc  = CriticalDensity()       ; %[kg/m^3]
    
    delta = rho / rhoc  ;
    tau   = Tc ./ T     ;
    
    cp = HeatCapacityIsobaricSingleR(delta,tau);
end