function [mayHaveSublimated,isBelowTriple,mayHaveFrozen] = mayHaveIceDeltaIND(delta,iND)

    %   Triple line values
    [delLt,delGt] = TriplePointDensitiesR();
    
    %   Ice-line masks
    canSublimate  = (delta >= +5.4715631599602330E-09) & (delta <= +1.5076377102923994E-05);
    canTriple     = (delta >=           delGt        ) & (delta <=           delLt        );
    canFreeze     = (delta >= +3.1049457143874224E+00) & (delta <= +3.3888854485186819E+00);

    %   Allocate masks
    deltaSize         = size(delta)     ;
    mayHaveSublimated = false(deltaSize);
    isBelowTriple     = false(deltaSize);
    mayHaveFrozen     = false(deltaSize);
    
    %   Evaluate
    iNDcompare                      = InternalEnergySublimationFitRND(delta(canSublimate))  ;
    mayHaveSublimated(canSublimate) = iNDcompare >= iND(canSublimate)                       ;
    iNDcompare                      = InternalEnergyTripleLineRND(delta(canTriple))         ;
    isBelowTriple(canTriple)        = iNDcompare >= iND(canTriple)                          ;
    iNDcompare                      = InternalEnergyFreezingFitRND(delta(canFreeze))        ;
    mayHaveFrozen(canFreeze)        = iNDcompare >= iND(canFreeze)                          ;


end