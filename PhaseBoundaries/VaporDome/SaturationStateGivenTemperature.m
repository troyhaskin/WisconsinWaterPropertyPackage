function [P,rhoL,rhoG] = SaturationStateGivenTemperature(T,varargin)
    
    if iscolumn(T)
        performTranspose = false();
    else
        T = T(:);
        performTranspose = true();
    end
    
    [~,rhoc,Tc] = Nondimensionalizers();
    tau = Tc./T;

    [P,delL,delG] = SaturationStateGivenTauRR(tau,varargin{:});
    rhoL = delL * rhoc;
    rhoG = delG * rhoc;
    
    if performTranspose
        P    = P'   ;
        rhoL = rhoL';
        rhoG = rhoG';
    end
    
end