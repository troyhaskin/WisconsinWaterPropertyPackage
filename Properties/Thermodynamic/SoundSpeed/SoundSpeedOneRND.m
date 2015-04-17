function wND = SoundSpeedOneRND(delta,tau)
    
    PhiR_d  = HelmholtzResidual_d (delta,tau)   ;
    PhiR_dd = HelmholtzResidual_dd(delta,tau)   ;
    PhiR_dt = HelmholtzResidual_dt(delta,tau)   ;
    Phi_tt  = Helmholtz_tt(delta,tau)           ;
    
    wND   = sqrt((1 + 2*delta.*PhiR_d + delta.^2.*PhiR_dd        -...
            (1 + delta.*PhiR_d - delta.*tau.*PhiR_dt).^2 ./...
            (tau.^2 .* Phi_tt)) ./ tau);
    
end