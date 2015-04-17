function Pnd_tau = PressureOneRND_tau(delta,tau)
    
    PhiR_d  = HelmholtzResidual_d (delta,tau);
    PhiR_dt = HelmholtzResidual_dt(delta,tau);
    
    Pnd_tau = (delta.*tau.*PhiR_dt - delta.*PhiR_d - 1) .* delta ./ tau.^2;

end