function s = EntropyOneR(delta,tau)
    
    sStar = DimensioningEntropy();
    sND   = EntropyOneRND(delta,tau);
    
    s = sND * sStar;
    
end