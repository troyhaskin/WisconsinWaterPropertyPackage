function [] = VerifyValues_HelmholtzResidual(rho,T)
    
    if (nargin == 0)    
    % IAPWS Test Values
        rho = 838.025   ; %[kg/m^3]
        T   = 500       ; %[K]
    end
   
    Tc    = CriticalTemperature()   ; %[K]
    rhoc  = CriticalDensity()       ; %[kg/m^3]
   
    
    delta = rho /  rhoc	;
    tau   = Tc  /   T   ;
    
   fprintf('\n');
   fprintf(' Sample Values for the Helmholtz Residual  functions \n');
   fprintf('=====================================================\n');
   fprintf('        State: {rho = %g, T = %g} \n',rho,T             );
   fprintf('-----------------------------------------------------\n');
   fprintf('PhiR    = %+17.8E\n',HelmholtzResidual   (delta,tau));
   fprintf('PhiR_d  = %+17.8E\n',HelmholtzResidual_d (delta,tau));
   fprintf('PhiR_dd = %+17.8E\n',HelmholtzResidual_dd(delta,tau));
   fprintf('PhiR_t  = %+17.8E\n',HelmholtzResidual_t (delta,tau));
   fprintf('PhiR_tt = %+17.8E\n',HelmholtzResidual_tt(delta,tau));
   fprintf('PhiR_dt = %+17.8E\n',HelmholtzResidual_dt(delta,tau));
   fprintf('=====================================================\n');
   fprintf('\n');
   
    
end