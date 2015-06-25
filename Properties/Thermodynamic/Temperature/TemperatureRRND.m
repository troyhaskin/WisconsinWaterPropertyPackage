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
        state2.isTwoPhi               = state.isTwoPhi  ;
        state                         = state2          ;
        
        %   Detect and contract out false-positives
        notTwoPhase             = xor(state.isTwoPhi,mayBeTwoPhase) ;
        notTwoPhase             = notTwoPhase (mayBeTwoPhase)       ;
        state.tau(notTwoPhase)  = []                                ;
        state.Pnd(notTwoPhase)  = []                                ;
        state.delL(notTwoPhase) = []                                ;
        state.delG(notTwoPhase) = []                                ;
        state.x(notTwoPhase)    = []                                ;
    end
    
end


