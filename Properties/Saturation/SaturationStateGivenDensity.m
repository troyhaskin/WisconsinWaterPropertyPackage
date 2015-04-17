function [P,T,rhoL,rhoG] = SaturationStateGivenDensity(rhoGiven,UniqueMask)
    
    
    if (nargin < 2)
        UniqueMask = [];
    end
    
    %   Non-dimensionalize
    rhoc     = CriticalDensity();
    delGiven = rhoGiven/rhoc;

    %   Call core function
    [Pnd,tau,delL,delG] = SaturationStateGivenDeltaRRND(delGiven,UniqueMask);
    
    %   Dimensionalize
    P    = Pnd * DimensioningPressure() ;
    T    = CriticalTemperature()./tau   ;
    rhoL = delL * rhoc                  ;
    rhoG = delG * rhoc                  ;
    
end



