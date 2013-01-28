function Deltabi_t = GetDeltabi_t(Delta,Theta,b)
    Deltabi_t   = -2*b * Theta .* Delta.^(b-1);
end
