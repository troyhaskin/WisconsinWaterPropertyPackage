function i_delta = InternalEnergyOneR_delta(delta,tau)
    
    istar     = DimensioningInternalEnergy();
    iND_delta = InternalEnergyOneRND_delta(delta,tau);
    
    i_delta = iND_delta * istar;
    
end