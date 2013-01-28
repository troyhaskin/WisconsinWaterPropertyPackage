function IntEnergy_tau = InternalEnergyOneR_tau(delta,tau)
    
    IntEnergyStar   = DimensioningInternalEnergy()          ;
    IntEnergyND_tau = InternalEnergyOneRND_tau(delta,tau)   ;
    
    IntEnergy_tau = IntEnergyND_tau * IntEnergyStar;
    
end