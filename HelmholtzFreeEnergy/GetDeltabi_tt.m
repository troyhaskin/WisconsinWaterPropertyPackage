function Deltabi_tt = GetDeltabi_tt(Delta,Theta,b)
    Deltabi_tt  = 2*b * Delta.^(b-1) + 4*b*(b-1) * Theta.^2 .* Delta.^(b-2);
end
