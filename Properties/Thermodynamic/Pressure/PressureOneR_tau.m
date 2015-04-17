function P_tau = PressureOneR_tau(delta,tau)
    
    Pstar   = DimensioningPressure()        ;
    Pnd_tau = PressureOneRND_tau(delta,tau) ;
    
    P_tau = Pnd_tau * Pstar;

end