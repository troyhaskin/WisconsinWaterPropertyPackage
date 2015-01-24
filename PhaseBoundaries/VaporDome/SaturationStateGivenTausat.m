function [P,delL,delG] = SaturationStateGivenTausat(tau,varargin)
    
    [P,delL,delG] = SaturationStateGivenTauRR(tau,varargin{:});
    
end