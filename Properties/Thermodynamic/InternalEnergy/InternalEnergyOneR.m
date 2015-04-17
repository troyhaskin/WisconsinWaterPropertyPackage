function IntEnergy = InternalEnergyOneR(delta,tau)

    IntEnergyStar = DimensioningInternalEnergy();
    IntEnergyND   = InternalEnergyOneRND(delta,tau);
    
    IntEnergy = IntEnergyND * IntEnergyStar;
    
end