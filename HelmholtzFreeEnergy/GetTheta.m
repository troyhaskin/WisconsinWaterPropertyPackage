function Theta = GetTheta(deltaMod,tau,A,betaInv)
    Theta = (1 - tau) + A .* deltaMod.^(0.5*betaInv);
end