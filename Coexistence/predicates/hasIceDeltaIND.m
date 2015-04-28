function [hasSublimated,isBelowTriple,hasFrozen] = hasIceDeltaIND(delta,iND)

    %   Perform approximate look-up
    [mayHaveSublimated,isBelowTriple,mayHaveFrozen] = mayHaveIceDeltaIND(delta,iND);

    %    Allocate
    hasSublimated = mayHaveSublimated;
    hasFrozen     = mayHaveFrozen;
    
    %   Call coexistence calculators
    iNDcompare                       = InternalEnergySublimationRND(delta(mayHaveSublimated))   ;
    hasSublimated(mayHaveSublimated) = iNDcompare >= iND(mayHaveSublimated)                     ;
    iNDcompare                       = InternalEnergyMeltRND(delta(mayHaveFrozen))              ;
    hasFrozen(mayHaveFrozen)         = iNDcompare >= iND(mayHaveFrozen)                         ;

end