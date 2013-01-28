function Helm = Helmholtz(delta,tau)
    Helm = HelmholtzIdealGas(delta,tau) + HelmholtzResidual(delta,tau);
end
