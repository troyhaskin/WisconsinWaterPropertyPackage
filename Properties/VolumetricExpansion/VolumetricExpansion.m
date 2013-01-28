function beta = VolumetricExpansion(rho,T,PhaseCheck,TwoPhiHandle)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 4) || isempty(TwoPhiHandle)
        TwoPhiHandle = 'void';
    end
    
    OnePhiHandle    = @(delta,tau,~) VolumetricExpansionOneR(delta,tau);
    
    beta = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiHandle,PhaseCheck);
    beta = RestoreShape(beta,GreatestProduct(SizeRho,SizeT));
    
end