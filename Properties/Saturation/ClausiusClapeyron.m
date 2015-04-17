function P_Tsat = ClausiusClapeyron(Psat,Tsat,rhol,rhog,varargin)
    
    [R,rhoc,Tc] = Nondimensionalizers();
    
    % Reduce the variables
    delL   = rhol  / rhoc       ;
    delG   = rhog  / rhoc       ;
    tauSat = Tc   ./ Tsat       ;
    Pnd    = Psat ./ (R*rhoc*Tc);
    
    % varargin is an optional input argument for the liquid and gas internal 
    % energies, respectively, if they have already been calculated.
    P_Tsat = ClausiusClapeyronR(Pnd,tauSat,delL,delG,varargin{:});
    
end