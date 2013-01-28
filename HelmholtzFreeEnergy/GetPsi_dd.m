function Psi_dd = GetPsi_dd(deltaMod,Psi,C)
    Psi_dd  = -2*C * Psi .* (2*C*deltaMod-1);
end
