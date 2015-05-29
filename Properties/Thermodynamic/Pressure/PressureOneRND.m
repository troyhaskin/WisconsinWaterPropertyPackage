function Pnd = PressureOneRND(delta,tau)
    
    PhiR_d = HelmholtzResidual_d(delta,tau);
    
    Pnd = (delta + delta.^2 .* PhiR_d)./ tau;
    
end