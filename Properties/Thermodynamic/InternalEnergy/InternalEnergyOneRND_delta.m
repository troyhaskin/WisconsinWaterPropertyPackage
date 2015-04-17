function iND_delta = InternalEnergyOneRND_delta(delta,tau)
    
    Phi_dt    = Helmholtz_dt(delta,tau);
    iND_delta = Phi_dt;
    
end