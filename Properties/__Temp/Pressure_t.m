function P = Pressure_t(rho,T,PhaseCheck)
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    R     = SpecificGasConstant()   ; %[J/kg/K]
    
    Rrho            = R .* rho;
    OnePhiHandle    = @(delta,tau,Mask) SinglePhasePressure_t(delta,tau,Mask,Rrho);            
	TwoPhiHandle	= @(~,~,~,~,~,TwoPhase) SmartMask(0,TwoPhase);
    
    P = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiHandle,PhaseCheck);
    
end

function Pone = SinglePhasePressure_t(delta,tau,Mask,Rrho)
    Pone = 1 +        delta .* HelmholtzResidual_d (delta,tau) + ...
               tau .* delta .* HelmholtzResidual_dt(delta,tau); % dimensionless
    Pone = Pone .* SmartMask(Rrho,Mask);                        % dimensional
end