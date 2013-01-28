function Psi_t = GetPsi_t(tau,Psi,D)
    Psi_t = -2*D * (tau-1) .* Psi;
end
