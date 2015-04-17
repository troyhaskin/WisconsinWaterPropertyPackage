function sND = EntropyOneRND(delta,tau)
    
    Phi   = Helmholtz  (delta,tau) ;
    Phi_t = Helmholtz_t(delta,tau) ;
    
    sND = tau.* Phi_t - Phi;
    
end