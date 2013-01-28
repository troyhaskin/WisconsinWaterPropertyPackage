function tau_T = tauRR_Temperature(tau)
    
    Tc    = CriticalTemperature();
    tau_T = -tau.^2 / Tc;
    
end