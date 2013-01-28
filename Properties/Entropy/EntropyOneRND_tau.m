function sND_tau = EntropyOneRND_tau(delta,tau)
    
    Phi_tt = Helmholtz_tt(delta,tau) ;
    
    sND_tau = tau.* Phi_tt;
    
end