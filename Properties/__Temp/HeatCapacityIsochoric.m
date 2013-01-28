function cv = HeatCapacityIsochoric(rho,T)
    
    Tc    = 647.096      ; %[K]
    rhoc  = 322          ; %[kg/m^3]
    R     = 0.46151805E3 ; %[J/kg/K]
    
    delta = rho./rhoc   ;
    tau   = Tc ./T      ;
    
    Phi0_tt = HelmholtzIdealGas_tt(delta,tau);
    PhiR_tt = HelmholtzResidual_tt(delta,tau);
    
    cv = -tau.^2 .* (Phi0_tt + PhiR_tt) * R;
end