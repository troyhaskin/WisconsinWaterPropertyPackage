function Pnd = PressureOneRND(delta,tau)
    
    PhiR_d = HelmholtzResidual_d(delta,tau);
    
    Pnd = (1 + delta .* PhiR_d) .* (delta ./ tau);
    
end