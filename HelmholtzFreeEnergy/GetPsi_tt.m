function Psi_tt = GetPsi_tt(tau,Psi,D)
    Psi_tt = 2*D * Psi .* (2*D*(tau-1).^2-1);
end
