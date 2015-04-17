function cpND = HeatCapacityIsobaricSingleRND(delta,tau)
    % Single phase property with reduced state specification
    % and non-dimensional output.
    
    Phi_tt  = Helmholtz_tt        (delta,tau);
    PhiR_d  = HelmholtzResidual_d (delta,tau);
    PhiR_dd = HelmholtzResidual_dd(delta,tau);
    PhiR_dt = HelmholtzResidual_dt(delta,tau);
    
    cvND = -tau.^2 .* Phi_tt;
    Rmod = (1 +   delta.*PhiR_d - delta.*tau.*PhiR_dt).^2 ./     ...
           (1 + 2*delta.*PhiR_d + delta.^2  .*PhiR_dd)           ;
    cpND = (cvND + Rmod);
end