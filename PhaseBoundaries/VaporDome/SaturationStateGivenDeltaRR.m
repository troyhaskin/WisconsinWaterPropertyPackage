function [P,tau,delL,delG] = SaturationStateGivenDeltaRR(delGiven,UniqueMask)
    
    if (nargin < 2)
        UniqueMask = [];
    end

    [Pnd,tau,delL,delG] = SaturationStateGivenDeltaRRND(delGiven,UniqueMask);
    P = Pnd * DimensioningPressure();

end







