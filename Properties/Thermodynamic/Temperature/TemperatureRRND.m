function [tau,state] = TemperatureRRND(delta,iND,tau0,PhaseCheck)
    
    %   Approximate phase check
    if PhaseCheck
        mayBeSaturated = mayBeSaturatedDeltaIND(delta,iND);
        mayHaveIce     = mayHaveIceDeltaIND(delta,iND)    ;
        
        onePhase      = not(mayBeSaturated) & not(mayHaveIce);
        mayBeTwoPhase = mayBeSaturated      & not(mayHaveIce); % false-positives handled in two-phase block
        
    else
        onePhase      = (delta == delta);
        mayBeTwoPhase = not(onePhase)   ;
    end
    
    %   Allocation
    n                   = numel(delta) ;
    tau(n,1)            = 0            ;
    state.isTwoPhi(n,1) = false        ;
    state.tau(n,1)      = 0            ;
    state.Pnd(n,1)      = 0            ;
    state.del(n,1)      = 0            ;
    state.delL(n,1)     = 0            ;
    state.delG(n,1)     = 0            ;
    state.x(n,1)        = -100         ;


    %   Single phase solve
    if any(onePhase)
        tau(onePhase) = TemperatureOneRRND(delta(onePhase),iND(onePhase),tau0(onePhase));
    end

    %   Two phase solve
    if any(mayBeTwoPhase)
        
        %   Calculate
        state2 = SaturationStateGivenMixedRhoIRRND(delta(mayBeTwoPhase),iND(mayBeTwoPhase),tau0(mayBeTwoPhase));
        
        %   Temperatures
        tau(mayBeTwoPhase) = state2.tau;
        
        %   Expand predicate to all states given to function
        state.isTwoPhi(mayBeTwoPhase) = state2.isTwoPhi ;
        state.tau(mayBeTwoPhase)      = state2.tau      ;
        state.Pnd(mayBeTwoPhase)      = state2.Pnd      ;
        state.del(mayBeTwoPhase)      = state2.delL     ;
        state.delL(mayBeTwoPhase)     = state2.delL     ;
        state.delG(mayBeTwoPhase)     = state2.delG     ;
        state.x(mayBeTwoPhase)        = state2.x        ;
        
    end
    
end


