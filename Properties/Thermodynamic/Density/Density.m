function rho = Density(P,T,varargin)

    delta = DensityRRND(P/DimensioningPressure(),CriticalTemperature()./T,varargin{:})  ;
    rho   = delta * CriticalDensity()                                       ;

end