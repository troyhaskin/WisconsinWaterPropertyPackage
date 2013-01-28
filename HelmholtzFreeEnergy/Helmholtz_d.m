function HelmDeriv = Helmholtz_d(delta,tau)
    HelmDeriv = HelmholtzIdealGas_d(delta,tau) + ...
                HelmholtzResidual_d(delta,tau);
end