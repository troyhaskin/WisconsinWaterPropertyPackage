function lamND = ThermalConductivityOneRND(delta,tau)
    
    lam0 = ThermalConductivityDiluteGasLimit           (tau);
    lam1 = ThermalConductivityFiniteDensity      (delta,tau);
    lam2 = ThermalConductivityCriticalEnhancement(delta,tau);
    
    lamND = lam0 .* lam1 + lam2;
    
end
