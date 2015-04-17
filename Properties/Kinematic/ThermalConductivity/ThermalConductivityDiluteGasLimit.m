function lam0 = ThermalConductivityDiluteGasLimit(tau)
    
    L        = ThermalConductivityDiluteGasLimitConstants();
    PowerSum = HornersMethod(tau,L);
    
    lam0 = 1 ./ (sqrt(tau) .* PowerSum);
    
end