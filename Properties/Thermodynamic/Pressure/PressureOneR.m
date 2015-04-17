function P = PressureOneR(delta,tau)
    
    Pstar = DimensioningPressure()      ;
    Pnd   = PressureOneRND(delta,tau)   ;
    
    P = Pnd * Pstar;
    
end