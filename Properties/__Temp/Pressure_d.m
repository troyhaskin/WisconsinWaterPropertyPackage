function P = Pressure_d(rho,T,PhaseCheck)
    
    if nargin < 3
        PhaseCheck = true;
    end
    
%     Tc    = CriticalTemperature()   ; %[K]
%     rhoc  = CriticalDensity()       ; %[kg/m^3]
    R     = SpecificGasConstant()   ; %[J/kg/K]
%     P     = Zeroer*0                ; % Allocate
    
    RT              = R .* T;
    OnePhiHandle    = @(delta,tau,Mask) SinglePhasePressure_d(delta,tau,Mask,RT);            
	TwoPhiHandle	= @(~,~,~,~,~,TwoPhase) SmartMask(0,TwoPhase);
    
    P = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiHandle,PhaseCheck);
    
end

function Pone = SinglePhasePressure_d(delta,tau,Mask,RT)
    Pone = 1 + delta    .* HelmholtzResidual_d (delta,tau) + ...
               delta.^2 .* HelmholtzResidual_dd(delta,tau); % dimensionless
    Pone = Pone .* SmartMask(RT,Mask);               % dimensional
end