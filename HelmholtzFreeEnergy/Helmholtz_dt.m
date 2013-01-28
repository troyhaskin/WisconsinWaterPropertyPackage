function HelmDeriv = Helmholtz_dt(delta,tau)
    HelmDeriv = HelmholtzIdealGas_dt(delta,tau) + ...
                HelmholtzResidual_dt(delta,tau);
end