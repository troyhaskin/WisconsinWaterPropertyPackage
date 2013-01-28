function Psi_d = GetPsi_d(delta,Psi,C)
    Psi_d  = -2*C * (delta-1) .* Psi;
end
