function [P,rhoL,rhoG] = SaturationStateGivenTemperature(T,varargin)
    
    [~,rhoc,Tc] = Nondimensionalizers();
    tau = Tc./T;

    [P,delL,delG] = SaturationStateGivenTauRRND(tau,varargin{:});
    rhoL = delL * rhoc;
    rhoG = delG * rhoc;
    
end