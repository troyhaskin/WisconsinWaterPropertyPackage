function delG_tau = SatGDensityRND_tau(tau,Psat,delL,delG,il,ig)
    
    % If the saturation state was not specified, calculate it.
    if (nargin < 2)
        [Psat,delL,delG] = SaturationStateGivenTausat(tau);
    end

    % If the internal energies were not specified, calculate them.
    if (nargin < 5)
        il = InternalEnergyOneR(delL,tau);
        ig = InternalEnergyOneR(delG,tau);
    end
    
    % Assign a generic del variable
    del = delG;
    
    % Dimensionless saturation Pressure and it's derivative w.r.t. tauSat
    Pstar   = DimensioningPressure();
    Pnd     = Psat / Pstar;
    Pnd_tau = ClausiusClapeyronRRND(Psat,tau,delL,delG,il,ig);
    
    % Helmholtz functions
    PhiR_d  = HelmholtzResidual_d (del,tau);
    PhiR_dt = HelmholtzResidual_dt(del,tau);
    PhiR_dd = HelmholtzResidual_dd(del,tau);
    
    % Derivative
    delG_tau = (del .* (Pnd + tau.*Pnd_tau - del.^2.*PhiR_dt)) ./ ...
               ( tau.*Pnd + del.^2.*(PhiR_d + del.*PhiR_dd) )  ;
    
end