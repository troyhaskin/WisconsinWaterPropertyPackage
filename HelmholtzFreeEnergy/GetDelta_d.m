function Delta_d = GetDelta_d(delta,deltaMod,Theta,A,B,a,betaInv)
    Part        = 2 .* betaInv .* A .* Theta .* deltaMod.^(0.5*betaInv - 1);
    Part        = Part + 2 * B .* a .* deltaMod.^(a-1);
    Delta_d     = (delta-1) .* Part;
end
