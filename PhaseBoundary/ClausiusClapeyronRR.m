function P_tauSat = ClausiusClapeyronRR(Psat,tauSat,delL,delG,varargin)
    
    Pstar = DimensioningPressure();
    
    % varargin is an optional input argument for the liquid and gas internal 
    % energies, respectively, if they have already been calculated.
    Pnd_tauSat = ClausiusClapeyronRRND(Psat,tauSat,delL,delG,varargin{:});

    % Chain rule
    P_tauSat = Pnd_tauSat * Pstar;
    
end