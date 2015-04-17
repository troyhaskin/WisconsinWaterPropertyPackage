function alpha = HomogeneousVoidFraction(rhoMix,rhoL,rhoG)

    alpha = (rhoMix - rhoL)./(rhoG - rhoL);
    
end