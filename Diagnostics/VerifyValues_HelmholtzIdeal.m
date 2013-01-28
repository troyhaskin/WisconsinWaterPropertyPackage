function [] = VerifyValues_HelmholtzIdeal()
   
    Tc    = 647.096      ; %[K]
    rhoc  = 322          ; %[kg/m^3]
    
    TTest   = 500       ; %[K]
    rhoTest = 838.025   ; %[kg/m^3]
    
    delta = rhoTest/  rhoc  ;
    tau   =   Tc   / TTest  ;
    
   fprintf('\n');
   fprintf(' Sample Values for the Helmholtz Ideal Gas functions \n');
   fprintf('=====================================================\n');
   fprintf('        State: {rho = %g, T = %g} \n',rhoTest,TTest     );
   fprintf('-----------------------------------------------------\n');
   fprintf('Phi0    = %+17.8E\n',HelmholtzIdealGas   (delta,tau));
   fprintf('Phi0_d  = %+17.8E\n',HelmholtzIdealGas_d (delta,tau));
   fprintf('Phi0_dd = %+17.8E\n',HelmholtzIdealGas_dd(delta,tau));
   fprintf('Phi0_t  = %+17.8E\n',HelmholtzIdealGas_t (delta,tau));
   fprintf('Phi0_tt = %+17.8E\n',HelmholtzIdealGas_tt(delta,tau));
   fprintf('Phi0_dt = %+17.8E\n',HelmholtzIdealGas_dt(delta,tau));
   fprintf('=====================================================\n');
   fprintf('\n');
   
    
end