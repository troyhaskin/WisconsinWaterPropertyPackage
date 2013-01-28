function tau_T = tau_TemperatureRR(tau)
    
    Tc    = CriticalTemperature();
    tau_T = -tau.^2 / Tc;
    
end