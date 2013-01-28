function i_rho = InternalEnergy_Density(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if nargin < 3
        PhaseCheck = true;
    end
       
    i_delta   = InternalEnergy_delta(rho,T,PhaseCheck);
    delta_rho = delta_Density();
    
    i_rho = i_delta .* delta_rho;
    i_rho = RestoreShape(i_rho,GreatestProduct(SizeRho,SizeT));
    
end