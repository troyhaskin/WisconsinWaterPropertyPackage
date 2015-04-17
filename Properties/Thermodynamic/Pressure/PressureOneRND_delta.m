function Pnd_delta = PressureOneRND_delta(delta,tau)
    
    PhiR_d  = HelmholtzResidual_d (delta,tau);
    PhiR_dd = HelmholtzResidual_dd(delta,tau);
    
    Pnd_delta = (1 + 2*delta.*PhiR_d + delta.^2 .* PhiR_dd) ./ tau;
    
end
