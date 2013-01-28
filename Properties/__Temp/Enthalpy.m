function h = Enthalpy(rho,T)
    
    Tc    = 647.096      ; %[K]
    rhoc  = 322          ; %[kg/m^3]
    R     = 0.46151805E3 ; %[J/kg/K]
    
    delta = rho./rhoc   ;
    tau   = Tc ./T      ;
    
    Phi0_t = HelmholtzIdealGas_t(delta,tau);
    PhiR_t = HelmholtzResidual_t(delta,tau);
    PhiR_d = HelmholtzResidual_d(delta,tau);
    
    h = (1 + tau .*(Phi0_t + PhiR_t) + delta.*PhiR_d) .* T * R;
end