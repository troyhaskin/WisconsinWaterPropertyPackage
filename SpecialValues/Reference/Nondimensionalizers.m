function [R,rhoc,Tc] = Nondimensionalizers()
   
    R    = SpecificGasConstant(); %[J/kg-K]
    rhoc = CriticalDensity    (); %[kg/m^3]
    Tc   = CriticalTemperature(); %[K]
    
end