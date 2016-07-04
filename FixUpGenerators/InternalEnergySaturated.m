function iND = InternalEnergySaturated(delta)
    
    [~,tau,~,~] = SaturationStateGivenDeltaRRND(delta(:));
    iND         = InternalEnergyOneRND(delta(:),tau).';
end