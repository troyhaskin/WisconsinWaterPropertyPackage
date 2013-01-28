function IntEnergyMix = InternalEnergyMixtureR(delMix,delL,delG,tauSat)
    
    % Quality weight
    x = QualityFromDensity(delMix,delL,delG);
    
    % Saturated internal energy values
    IntEnergyMixL = InternalEnergyOneR(delL,tauSat);
	IntEnergyMixG = InternalEnergyOneR(delG,tauSat);
    
    % Mixture internal energy
	IntEnergyMix = IntEnergyMixL + x .* (IntEnergyMixG - IntEnergyMixL);

end