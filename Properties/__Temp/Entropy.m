function s = Entropy(rho,T)
    
    Tc    = 647.096      ; %[K]
    rhoc  = 322          ; %[kg/m^3]
    R     = 0.46151805E3 ; %[J/kg/K]
    
    delta = rho/rhoc;
    tau   = Tc./T;
    
    Phi0   = HelmholtzIdealGas  (delta,tau) ;
    PhiR   = HelmholtzResidual  (delta,tau) ;
    Phi0_t = HelmholtzIdealGas_t(delta,tau) ;
    PhiR_t = HelmholtzResidual_t(delta,tau) ;
    s      = (tau .* (Phi0_t + PhiR_t) - Phi0 - PhiR) * R;
end