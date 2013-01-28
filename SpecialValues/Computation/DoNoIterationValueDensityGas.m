function rhoG = DoNoIterationValueDensityGas()
    %   Any density above this value will use the curve fits and not 
    %   iteration for the saturation properties due to an inability to 
    %   get the system to robustly converge above the value;
    delDoNo = DoNoIterationValueReducedDensityGas() ;
    rhoc    = CriticalDensity()                     ;
    rhoG    = delDoNo * rhoc;
end