function Pnd_Tsat = ClausiusClapeyronRND(Psat,tauSat,delL,delG,varargin)
    
    % varargin is an optional input argument for the liquid and gas internal 
    % energies, respectively, if they have already been calculated.
    Pnd_tauSat = ClausiusClapeyronRRND(Psat,tauSat,delL,delG,varargin{:});
    T_tauSat   = Temperature_tau(tauSat);
    
    % Chain Rule
    Pnd_Tsat = Pnd_tauSat ./ T_tauSat;
    
end