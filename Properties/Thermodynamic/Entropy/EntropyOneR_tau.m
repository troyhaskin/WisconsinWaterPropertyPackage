function s_tau = EntropyOneR_tau(delta,tau)
    
    sStar   = DimensioningEntropy();
    sND_tau = EntropyOneRND_tau(delta,tau);
    
    s_tau = sND_tau * sStar;
    
end