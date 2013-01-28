function Deltabi_d = GetDeltabi_d(delta,deltaMod,Delta,Theta,A,B,a,b,betaInv)
    Delta_d     = GetDelta_d(delta,deltaMod,Theta,A,B,a,betaInv);
    Deltabi_d   = b * Delta.^(b-1).*Delta_d;
end
