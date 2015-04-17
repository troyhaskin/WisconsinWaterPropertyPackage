function P_Tsat = ClausiusClapeyronR(Pnd,tauSat,delL,delG,varargin)
    
    Pstar    = DimensioningPressure();
    Pnd_Tsat = ClausiusClapeyronRND(Pnd,tauSat,delL,delG,varargin{:});

    P_Tsat = Pnd_Tsat * Pstar;

end



