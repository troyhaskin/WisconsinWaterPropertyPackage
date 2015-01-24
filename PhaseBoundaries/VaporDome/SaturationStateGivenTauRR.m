function [P,delL,delG] = SaturationStateGivenTauRR(tau,varargin)

    [Pnd,delL,delG] = SaturationStateGivenTauRRND(tau,varargin{:});
    P = Pnd * DimensioningPressure();
    
end