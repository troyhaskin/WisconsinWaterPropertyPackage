function hND = EnthalpyOneRND(delta,tau)
    
    Phi_t  = Helmholtz_t        (delta,tau);
    PhiR_d = HelmholtzResidual_d(delta,tau);
    
    hND = (1+ tau .* Phi_t + delta .* PhiR_d) ./ tau;
    
end