function P_delta = PressureOneR_delta(delta,tau)
    
    Pstar     = DimensioningPressure(); %[Pa]    
    Pnd_delta = PressureOneRND_delta(delta,tau);
    
    P_delta = Pnd_delta * Pstar;
end