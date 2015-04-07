function iNDt = InternalEnergyTripleLineRND(delta)
    
    [delLt,delGt] = TriplePointDensitiesR();
    
    iNDt = delta*0;
    
    isMixed = (delta <= delLt) && (delta >= delGt);
    TriplePointInternalEnergiesND
    
    
    
end