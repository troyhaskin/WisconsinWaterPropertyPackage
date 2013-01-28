function dPddelta = Pressure_d(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3)
        PhaseCheck = true;
    end
    
    R     = SpecificGasConstant()   ; %[J/kg/K]
    rhoc  = CriticalDensity()       ; %[kg/m^3]
    
    RT              = rhoc * R .* T;
    OnePhiHandle    = @(delta,tau,Mask) SinglePhasePressure_d(delta,tau,Mask,RT);            
	TwoPhiHandle	= @(~,~,~,~,~,TwoPhase) SmartMask(0,TwoPhase);
    
    dPddelta = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiHandle,PhaseCheck);
    dPddelta = RestoreShape(dPddelta,GreatestProduct(SizeRho,SizeT));
end

function Pone = SinglePhasePressure_d(delta,tau,Mask,rhocRT)
    Pone = 1 + delta    .* HelmholtzResidual_d (delta,tau) + ...
               delta.^2 .* HelmholtzResidual_dd(delta,tau); % dimensionless
    Pone = Pone .* SmartMask(rhocRT,Mask);                  % dimensional
end