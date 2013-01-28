function IntEnergyStar = DimensioningInternalEnergy()

    [R,~,Tc]      = Nondimensionalizers();
    IntEnergyStar = R * Tc;

end