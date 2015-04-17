function s_tau = Entropy_tau(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3)
        PhaseCheck = true;
    end
    
    OnePhiHandle    = @(delta,tau,~) EntropyOneR_tau(delta,tau);
    
    s_tau = GenericTwoPhiProperty(rho,T,OnePhiHandle,'quality',PhaseCheck);
    s_tau = RestoreShape(s_tau,GreatestProduct(SizeRho,SizeT));
    
end