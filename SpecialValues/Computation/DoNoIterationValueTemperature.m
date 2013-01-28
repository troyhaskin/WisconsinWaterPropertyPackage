function T = DoNoIterationValueTemperature()
    %   Any temperature above this value will use the curve fits and not 
    %   iteration for the saturation properties due to an inability to 
    %   get the system to robustly converge above the value;
    tauDoNo = DoNoIterationValueReducedTemperature();
    Tc      = CriticalTemperature()                 ;
    
    T       = Tc / tauDoNo;

end