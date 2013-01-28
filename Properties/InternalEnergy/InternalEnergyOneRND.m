function IntEnergyND   = InternalEnergyOneRND(delta,tau)
    
    Phi_t = Helmholtz_t(delta,tau);
    
    IntEnergyND = Phi_t;
    
end