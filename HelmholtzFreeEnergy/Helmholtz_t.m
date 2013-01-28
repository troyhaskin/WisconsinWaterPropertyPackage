function HelmDeriv = Helmholtz_t(delta,tau)
    HelmDeriv = HelmholtzIdealGas_t(delta,tau) + ...
                HelmholtzResidual_t(delta,tau);
end