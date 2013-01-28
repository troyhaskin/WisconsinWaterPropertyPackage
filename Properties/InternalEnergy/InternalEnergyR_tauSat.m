function IntEnergy_tauSat = InternalEnergyR_tauSat(tau,delMix,Psat,delL,delG,il,ig)
    
    % If the saturation state is not passed, calculate it.
    if (nargin < 3)
        [Psat,delL,delG] = SaturationStateGivenTausat(tau);
    end
    
    % If internal energies are not passed, calculate them.
    if (nargin < 6)
        il = InternalEnergyOneR(delL,tau);
        ig = InternalEnergyOneR(delG,tau); 
    end
    
    % Calculate the quality
	x = QualityFromDensity(delMix,delL,delG);
    
    % Calculate the saturated density's derivative w.r.t. tau
    delL_tau = SatLDensityRND_tau(tau,Psat,delL,delG,il,ig);
    delG_tau = SatGDensityRND_tau(tau,Psat,delL,delG,il,ig);
    
    % Calculate the quality's derivative w.r.t. tau
    x_tau = QualityR_tau(tau,delMix,Psat,delL,delG,delL_tau,delG_tau,il,ig);

    % Internal Energy derivatives
    il_tau = InternalEnergyOneR_tau  (delL,tau) ;
    ig_tau = InternalEnergyOneR_tau  (delG,tau) ;
    il_del = InternalEnergyOneR_delta(delL,tau) ;
    ig_del = InternalEnergyOneR_delta(delG,tau) ;
    
    % Induced tau derivatives:
    %       In two-phase, the deltas are NOT independent of the tau; therefore, 
    %       the internal energy's delta inputs induce tau derivatives.
    il_del_tau = il_del.*delL_tau;
    ig_del_tau = ig_del.*delG_tau;
    
    %   Internal energy derivatives for tau.
    iMix_del = il_del_tau + x.*(ig_del_tau - il_del_tau);
    iMix_tau = il_tau     + x.*(ig_tau     - il_tau    );
    
    % Total mixture internal energy's derivative w.r.t. tau
    IntEnergy_tauSat = x_tau .* (ig - il) + iMix_tau + iMix_del;

end




