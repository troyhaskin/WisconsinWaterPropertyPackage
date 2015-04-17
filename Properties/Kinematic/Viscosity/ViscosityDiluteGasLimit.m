function mu0 = ViscosityDiluteGasLimit(tau)
    
    H        = ViscosityDiluteGasLimitConstants();
    PowerSum = HornersMethod(tau,H);
    
    mu0 = 100 ./ (sqrt(tau) .* PowerSum);
    
end