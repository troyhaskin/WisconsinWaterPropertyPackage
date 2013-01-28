function Psi_dt = GetPsi_dt(delta,tau,Psi,C,D)
    Psi_dt =4*C*D * (delta-1) .* (tau-1) .* Psi;
end
