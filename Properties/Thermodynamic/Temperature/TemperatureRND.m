function [T,state] = TemperatureRND(delta,iND,Tguess,PhaseCheck)

    Tc          = CriticalTemperature()                     ;
    tau         = Tc ./ Tguess                              ;
    [tau,state] = TemperatureRRND(delta,iND,tau,PhaseCheck) ;
    T           = Tc ./ tau                                 ;

end