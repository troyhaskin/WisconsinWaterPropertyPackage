function IntEnergyND_tau = InternalEnergyOneRND_tau(delta,tau)
    
    Phi_tt = Helmholtz_tt(delta,tau);
    
    IntEnergyND_tau = Phi_tt;
    
end