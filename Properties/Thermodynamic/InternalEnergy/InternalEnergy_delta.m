function i_delta = InternalEnergy_delta(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if nargin < 3
        PhaseCheck = true;
    end
       
    OnePhiHandle = @(delta,tau,~) InternalEnergyOneR_delta(delta,tau);
%     TwoPhiHandle = @(Psat,delL,delG,del,tau,~) InternalEnergyTwoR_tau(Psat,delL,delG,del,tau);
    
    i_delta = GenericTwoPhiProperty(rho,T,OnePhiHandle,'quality',PhaseCheck);
    i_delta = RestoreShape(i_delta,GreatestProduct(SizeRho,SizeT));
    
end