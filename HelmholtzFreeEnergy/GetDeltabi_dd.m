function Deltabi_dd = GetDeltabi_dd(delta,deltaMod,Delta,Theta,A,B,a,b,betaInv)
    Delta_d     = GetDelta_d (delta,deltaMod,Theta,A,B,a,betaInv);
    Delta_dd    = GetDelta_dd(delta,deltaMod,Theta,A,B,a,betaInv);
    Deltabi_dd  = b * (Delta.^(b-1).*Delta_dd + (b-1)*Delta.^(b-2).*Delta_d.^2) ;
end
