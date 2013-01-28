function Deltabi_dt = GetDeltabi_dt(delta,deltaMod,Delta,Theta,A,B,a,b,betaInv)
    Delta_d     = GetDelta_d(delta,deltaMod,Theta,A,B,a,betaInv);
    Part        = -2*A*b*betaInv * Delta.^(b-1) .* (delta-1) .* deltaMod.^(0.5*betaInv-1);
    Deltabi_dt  = Part - 2*b*(b-1) * Theta .* Delta.^(b-2) .* Delta_d;
end
