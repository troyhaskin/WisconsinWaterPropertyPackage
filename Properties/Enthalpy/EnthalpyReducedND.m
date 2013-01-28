function h = EnthalpyReducedND(delta,tau,varargin)
    
    Phi0_t = HelmholtzIdealGas_t(delta,tau);
    PhiR_t = HelmholtzResidual_t(delta,tau);
    PhiR_d = HelmholtzResidual_d(delta,tau);
    
    h = (1 + tau .*(Phi0_t + PhiR_t) + delta.*PhiR_d);
end