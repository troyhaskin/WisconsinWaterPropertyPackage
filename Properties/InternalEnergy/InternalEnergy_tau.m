function i_tau = InternalEnergy_tau(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if nargin < 3
        PhaseCheck = true;
    end
       
    OnePhiHandle = @(delta,tau,~) InternalEnergyOneR_tau(delta,tau);
    TwoPhiHandle = @(Psat,delL,delG,del,tau,~) InternalEnergyTwoR_tau(Psat,delL,delG,del,tau);
    
    i_tau = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiHandle,PhaseCheck);
    i_tau = RestoreShape(i_tau,GreatestProduct(SizeRho,SizeT));
    
end

function iMix_tau = InternalEnergyTwoR_tau(Psat,delL,delG,del,tau)
    
    il     = InternalEnergyOneR    (delL,tau);
    ig     = InternalEnergyOneR    (delG,tau);
    il_tau = InternalEnergyOneR_tau(delL,tau);
    ig_tau = InternalEnergyOneR_tau(delG,tau);
    
    x     = QualityFromDensity(del,delL,delG);
    x_tau = QualityR_tauSat(tau,del,Psat,delL,delG,il,ig);
    
    iMix_tau = il_tau + x.*(ig_tau-il_tau) + (ig-il).*x_tau;
end
