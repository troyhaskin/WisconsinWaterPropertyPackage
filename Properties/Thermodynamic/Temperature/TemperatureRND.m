function T = TemperatureRND(delta,iND,Tguess,PhaseCheck)
    
    Tc  = CriticalTemperature()                     ;
    tau = Tc ./ Tguess                              ;
    tau = TemperatureRRND(delta,iND,tau,PhaseCheck) ;
    T   = Tc ./ tau                                 ;
    
end