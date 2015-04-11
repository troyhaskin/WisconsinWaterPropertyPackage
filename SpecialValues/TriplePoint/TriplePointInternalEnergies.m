function [iLt,iGt] = TriplePointInternalEnergies()
    
    [iLtND,iGtND] = TriplePointInternalEnergiesND();
    iLt = iLtND * DimensioningInternalEnergy();
    iGt = iGtND * DimensioningInternalEnergy();

end