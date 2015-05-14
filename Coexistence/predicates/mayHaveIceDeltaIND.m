function varargout = mayHaveIceDeltaIND(delta,iND)

    %   Triple line values
    [delLt,delGt] = TriplePointDensitiesR();
    
    %   Ice-line masks
    canSublimate  = (delta >= +5.4715631599602330E-09) & (delta <= +1.5076377102923994E-05);
    canTriple     = (delta >=           delGt        ) & (delta <=           delLt        );
    canFreeze     = (delta >= +3.1049457143874224E+00) & (delta <= +3.3888854485186819E+00);

    %   Allocate masks
    deltaSize         = size(delta)         ;
    mayHaveSublimated = false(deltaSize)    ;
    isBelowTriple     = mayHaveSublimated   ;
    mayHaveFrozen     = mayHaveSublimated   ;
    
    %   Evaluate
    iNDcompare                      = InternalEnergySublimationFitRND(delta(canSublimate))  ;
    mayHaveSublimated(canSublimate) = iNDcompare >= iND(canSublimate)                       ;
    iNDcompare                      = InternalEnergyTripleLineRND(delta(canTriple))         ;
    isBelowTriple(canTriple)        = iNDcompare >= iND(canTriple)                          ;
    iNDcompare                      = InternalEnergyFreezingFitRND(delta(canFreeze))        ;
    mayHaveFrozen(canFreeze)        = iNDcompare >= iND(canFreeze)                          ;

    %   Main predicate
    mayHaveIce = mayHaveSublimated | isBelowTriple | mayHaveFrozen;
    
    %   Variable output
    switch(nargout)
        case(1)
            varargout{1} = mayHaveIce;

        case(2)
            varargout = {mayHaveIce , isBelowTriple} ;

        case(3)
            varargout = {mayHaveSublimated , isBelowTriple , mayHaveFrozen} ;

        case(4)
            varargout = {mayHaveIce , mayHaveSublimated , isBelowTriple , mayHaveFrozen} ;

        otherwise
            error('WWPP:mayHaveIceDeltaIND:WrongOutputCount',...
                'Only 1, 2, 3, or 4 outputs are allowed.');
    end

end