function alpha = VoidFractionFromDensity(rhoMix,rhoL,rhoG,SlipRatio)
    
    x   = QualityFromDensity(rhoMix,rhoL,rhoG)  ;
    Psi = (rhoG./rhoL) .* SlipRatio             ;
    
    alpha = x./(x + (1-x) .* Psi);
    
end