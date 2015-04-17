function cvND = HeatCapacityIsochoricOneRND(delta,tau)
    
    Phi_tt = Helmholtz_tt(delta,tau);
    
    cvND = -tau.^2 .* Phi_tt;
    
end