function HelmDeriv = Helmholtz_dd(delta,tau)
    HelmDeriv = HelmholtzIdealGas_dd(delta,tau) + ...
                HelmholtzResidual_dd(delta,tau);
end