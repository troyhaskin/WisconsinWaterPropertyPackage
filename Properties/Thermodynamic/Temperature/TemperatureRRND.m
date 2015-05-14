function tau = TemperatureRRND(delta,iND,tau0,PhaseCheck)
    
    %   Approximate phase check
    if PhaseCheck
        mayBeSaturated = mayBeSaturatedDeltaIND(delta,iND);
        mayHaveIce     = mayHaveIceDeltaIND(delta,iND)    ;
        
        onePhase = not(mayBeSaturated) & not(mayHaveIce);
        twoPhase = mayBeSaturated      & not(mayHaveIce); % false-positives detected in two-phase solver
        
    else
        onePhase = (delta == delta) ;
        twoPhase = not(onePhase)    ;
    end
    
    %   Allocation
    tau = delta * 0;
    
    %   Single phase solve
    if any(onePhase)
        tau(onePhase) = TemperatureOneRRND(delta,iND,tau0);
    end

    %   Single phase solve
    if any(twoPhase)
        [~,tau(twoPhase),~,~,~] = SaturationStateGivenMixedRhoIRRND(delta,iND,tau0);
    end
    
end


