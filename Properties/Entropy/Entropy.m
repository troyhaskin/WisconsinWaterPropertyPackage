function s = Entropy(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3)
        PhaseCheck = true;
    end
    
    OnePhiHandle    = @(delta,tau,~) EntropyOneR(delta,tau);
    
    s = GenericTwoPhiProperty(rho,T,OnePhiHandle,'void',PhaseCheck);
    s = RestoreShape(s,GreatestProduct(SizeRho,SizeT));
    
end