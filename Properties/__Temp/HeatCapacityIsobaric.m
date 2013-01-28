function cp = HeatCapacityIsobaric(rho,T)
    
    Tc    = 647.096      ; %[K]
    rhoc  = 322          ; %[kg/m^3]
    R     = 0.46151805E3 ; %[J/kg/K]
    
    delta = rho./rhoc   ;
    tau   = Tc ./T      ;
    
    Phi0_tt = HelmholtzIdealGas_tt(delta,tau);
    PhiR_tt = HelmholtzResidual_tt(delta,tau);
    PhiR_d  = HelmholtzResidual_d (delta,tau);
    PhiR_dd = HelmholtzResidual_dd(delta,tau);
    PhiR_dt = HelmholtzResidual_dt(delta,tau);
    
    cv  = -tau.^2 .* (Phi0_tt + PhiR_tt);
    FTR = (1 +   delta.*PhiR_d - delta.*tau.*PhiR_dt).^2 ./     ...
          (1 + 2*delta.*PhiR_d + delta.^2  .*PhiR_dd)           ;
    cp  = (cv + FTR) * R;
end