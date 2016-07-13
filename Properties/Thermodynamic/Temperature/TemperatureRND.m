function [T,state] = TemperatureRND(delta,iND,tauGuess,PhaseCheck)

    Tc          = CriticalTemperature()                         ;
    tau         = Tc ./ Tguess                                  ;
    [tau,state] = TemperatureRRND(delta,iND,tauGuess,PhaseCheck);
    T           = Tc ./ tau                                     ;

end