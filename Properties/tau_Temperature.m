function tau_T = tau_Temperature(T)
    
    Tc    = CriticalTemperature()   ;
    tau_T = -Tc ./ T.^2             ;
    
end