function Psi = GetPsi(deltaMod,tau,C,D)
    Psi = exp(-C*deltaMod-D*(tau-1).^2);
end
