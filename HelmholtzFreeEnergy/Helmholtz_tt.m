function HelmDeriv = Helmholtz_tt(delta,tau)
    HelmDeriv = HelmholtzIdealGas_tt(delta,tau) + ...
                HelmholtzResidual_tt(delta,tau);
end