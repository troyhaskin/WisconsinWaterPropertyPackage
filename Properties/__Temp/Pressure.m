function P = Pressure(rho,T,PhaseCheck)
    
    [rho,T,SizeRho,SizeT] = PropertyReshape(rho,T);
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    R     = SpecificGasConstant()   ; %[J/kg/K]
    
    rhoRT           = rho .* R .* T;
    OnePhiHandle    = @(delta,tau,Mask) SinglePhasePressure(delta,tau,Mask,rhoRT);            
	TwoPhiHandle	= @(~,~,Psat,~,~,TwoPhase) SmartMask(Psat,TwoPhase);
    
    P = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiHandle,PhaseCheck);
    
    P = PropertyReshape(P,SizeRho,SizeT);
    
end

function Pone = SinglePhasePressure(delta,tau,Mask,rhoRT)
    Pone = 1 + delta .* HelmholtzResidual_d(delta,tau); % dimensionless
    Pone = Pone .* SmartMask(rhoRT,Mask);               % dimensional
end