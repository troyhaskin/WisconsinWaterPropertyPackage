function P_Tsat = ClausiusClapeyronR(Psat,tauSat,delL,delG,varargin)
    
    Pstar    = DimensioningPressure();
    Pnd_Tsat = ClausiusClapeyronRND(Psat,tauSat,delL,delG,varargin);

    P_Tsat = Pnd_Tsat * Pstar;

end



