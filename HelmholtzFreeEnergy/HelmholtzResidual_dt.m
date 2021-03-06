function HelmDeriv = HelmholtzResidual_dt(delta,tau)
    
    [c,d,t,n,alpha,beta,gamma,epsilon,A,B,C,D,a,b] = Coefficients_HelmholtzResidual();
    
    betaInv   = 1./beta         ;
	deltaMod  = (delta - 1).^2 + eps(delta)	;
    HelmDeriv = 0*delta               ;
    SumErr    = 0*delta               ; % Summation correction (Kahan summation)

    
    for k = 1:7
        Part      = n(k)*d(k)*t(k) * delta.^(d(k)-1) .* tau.^(t(k)-1)   ;
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
    
    for k = 8:51
        delta2c   = delta.^(c(k));
        Part      = (d(k)-c(k)*delta2c) .* exp(-delta2c)               	;
        Part      = n(k)*t(k) * delta.^(d(k)-1) .* tau.^(t(k)-1) .* Part;
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
    
    for k = 52:54
        m       = k - 51                                                        ;
        Arg     = -alpha(m)*(delta-epsilon(m)).^2 - beta(m)*(tau-gamma(m)).^2   ;
        Part    =          t(k)./tau   - 2*beta(m) *(tau  -gamma(m)   )         ;
        Part    = Part .* (d(k)./delta - 2*alpha(m)*(delta-epsilon(m)))         ;
        Part    = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(Arg) .* Part       ;
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
    
    for k = 55:56
        m        = k - 51                                       ;
        p        = k - 54                                       ;

        Theta       = GetTheta     (deltaMod,tau,A(p),betaInv(m))  ;
        Delta       = GetDelta     (deltaMod,Theta,B(p),a(p))   ;
        Deltabi_d   = GetDeltabi_d (delta,deltaMod,Delta,Theta,A(p),B(p),a(p),b(p),betaInv(m));
        Deltabi_t   = GetDeltabi_t (Delta,Theta,b(p))           ;
        Deltabi_dt  = GetDeltabi_dt(delta,deltaMod,Delta,Theta,A(p),B(p),a(p),b(p),betaInv(m));
        Psi         = GetPsi       (deltaMod,tau,C(p),D(p))     ;
        Psi_d       = GetPsi_d     (delta,Psi,C(p))             ;
        Psi_t       = GetPsi_t     (tau  ,Psi,D(p))             ;
        Psi_dt      = GetPsi_dt    (delta,tau,Psi,C(p),D(p))    ;
        
        Part        = Delta.^(b(p)) .* (Psi_t + delta .* Psi_dt);
        Part        = Part + delta .* Deltabi_d .* Psi_t        ;
        Part        = Part + Deltabi_t .*(Psi + delta .* Psi_d) ;
        Part        = n(k)*(Part + Deltabi_dt .* delta .* Psi)  ; 
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
    
end

% 
% function [c,d,t,n] = Coefficients_Group1()
%     
%     c(1)  = 0.0;
%     c(2)  = 0.0;
%     c(3)  = 0.0;
%     c(4)  = 0.0;
%     c(5)  = 0.0;
%     c(6)  = 0.0;
%     c(7)  = 0.0;
%     c(8)  = 1.0;
%     c(9)  = 1.0;
%     c(10) = 1.0;
%     c(11) = 1.0;
%     c(12) = 1.0;
%     c(13) = 1.0;
%     c(14) = 1.0;
%     c(15) = 1.0;
%     c(16) = 1.0;
%     c(17) = 1.0;
%     c(18) = 1.0;
%     c(19) = 1.0;
%     c(20) = 1.0;
%     c(21) = 1.0;
%     c(22) = 1.0;
%     c(23) = 2.0;
%     c(24) = 2.0;
%     c(25) = 2.0;
%     c(26) = 2.0;
%     c(27) = 2.0;
%     c(28) = 2.0;
%     c(29) = 2.0;
%     c(30) = 2.0;
%     c(31) = 2.0;
%     c(32) = 2.0;
%     c(33) = 2.0;
%     c(34) = 2.0;
%     c(35) = 2.0;
%     c(36) = 2.0;
%     c(37) = 2.0;
%     c(38) = 2.0;
%     c(39) = 2.0;
%     c(40) = 2.0;
%     c(41) = 2.0;
%     c(42) = 2.0;
%     c(43) = 3.0;
%     c(44) = 3.0;
%     c(45) = 3.0;
%     c(46) = 3.0;
%     c(47) = 4.0;
%     c(48) = 6.0;
%     c(49) = 6.0;
%     c(50) = 6.0;
%     c(51) = 6.0;
%     
%     d(1)  =   1.0;
%     d(2)  =   1.0;
%     d(3)  =   1.0;
%     d(4)  =   2.0;
%     d(5)  =   2.0;
%     d(6)  =   3.0;
%     d(7)  =   4.0;
%     d(8)  =   1.0;
%     d(9)  =   1.0;
%     d(10) =   1.0;
%     d(11) =   2.0;
%     d(12) =   2.0;
%     d(13) =   3.0;
%     d(14) =   4.0;
%     d(15) =   4.0;
%     d(16) =   5.0;
%     d(17) =   7.0;
%     d(18) =   9.0;
%     d(19) =  10.0;
%     d(20) =  11.0;
%     d(21) =  13.0;
%     d(22) =  15.0;
%     d(23) =   1.0;
%     d(24) =   2.0;
%     d(25) =   2.0;
%     d(26) =   2.0;
%     d(27) =   3.0;
%     d(28) =   4.0;
%     d(29) =   4.0;
%     d(30) =   4.0;
%     d(31) =   5.0;
%     d(32) =   6.0;
%     d(33) =   6.0;
%     d(34) =   7.0;
%     d(35) =   9.0;
%     d(36) =   9.0;
%     d(37) =   9.0;
%     d(38) =   9.0;
%     d(39) =   9.0;
%     d(40) =  10.0;
%     d(41) =  10.0;
%     d(42) =  12.0;
%     d(43) =   3.0;
%     d(44) =   4.0;
%     d(45) =   4.0;
%     d(46) =   5.0;
%     d(47) =  14.0;
%     d(48) =   3.0;
%     d(49) =   6.0;
%     d(50) =   6.0;
%     d(51) =   6.0;
%     d(52) =   3.0;
%     d(53) =   3.0;
%     d(54) =   3.0;
%     
%     t(1)  =   -0.5   ;
%     t(2)  =    0.875 ;
%     t(3)  =    1.0   ;
%     t(4)  =    0.5   ;
%     t(5)  =    0.75  ;
%     t(6)  =    0.375 ;
%     t(7)  =    1.0   ;
%     t(8)  =    4.0   ;
%     t(9)  =    6.0   ;
%     t(10) =   12.0   ;
%     t(11) =    1.0   ;
%     t(12) =    5.0   ;
%     t(13) =    4.0   ;
%     t(14) =    2.0   ;
%     t(15) =   13.0   ;
%     t(16) =    9.0   ;
%     t(17) =    3.0   ;
%     t(18) =    4.0   ;
%     t(19) =   11.0   ;
%     t(20) =    4.0   ;
%     t(21) =   13.0   ;
%     t(22) =    1.0   ;
%     t(23) =    7.0   ;
%     t(24) =    1.0   ;
%     t(25) =    9.0   ;
%     t(26) =   10.0   ;
%     t(27) =   10.0   ;
%     t(28) =    3.0   ;
%     t(29) =    7.0   ;
%     t(30) =   10.0   ;
%     t(31) =   10.0   ;
%     t(32) =    6.0   ;
%     t(33) =   10.0   ;
%     t(34) =   10.0   ;
%     t(35) =    1.0   ;
%     t(36) =    2.0   ;
%     t(37) =    3.0   ;
%     t(38) =    4.0   ;
%     t(39) =    8.0   ;
%     t(40) =    6.0   ;
%     t(41) =    9.0   ;
%     t(42) =    8.0   ;
%     t(43) =   16.0   ;
%     t(44) =   22.0   ;
%     t(45) =   23.0   ;
%     t(46) =   23.0   ;
%     t(47) =   10.0   ;
%     t(48) =   50.0   ;
%     t(49) =   44.0   ;
%     t(50) =   46.0   ;
%     t(51) =   50.0   ;
%     t(52) =    0.0   ;
%     t(53) =    1.0   ;
%     t(54) =    4.0   ;
%     
%     n(1)  =  0.12533547935523E-1 ;
%     n(2)  =  0.78957634722828E+1 ;
%     n(3)  = -0.87803203303561E+1 ;
%     n(4)  =  0.31802509345418    ;
%     n(5)  = -0.26145533859358    ;
%     n(6)  = -0.78199751687981E-2 ;
%     n(7)  =  0.88089493102134E-2 ;
%     n(8)  = -0.66856572307965    ;
%     n(9)  =  0.20433810950965    ;
%     n(10) = -0.66212605039687E-4 ;
%     n(11) = -0.19232721156002    ;
%     n(12) = -0.25709043003438    ;
%     n(13) =  0.16074868486251    ;
%     n(14) = -0.40092828925807E-1 ;
%     n(15) =  0.39343422603254E-6 ;
%     n(16) = -0.75941377088144E-5 ;
%     n(17) =  0.56250979351888E-3 ;
%     n(18) = -0.15608652257135E-4 ;
%     n(19) =  0.11537996422951E-8 ;
%     n(20) =  0.36582165144204E-6 ;
%     n(21) = -0.13251180074668E-11;
%     n(22) = -0.62639586912454E-9 ;
%     n(23) = -0.10793600908932    ;
%     n(24) =  0.17611491008752E-1 ;
%     n(25) =  0.22132295167546    ;
%     n(26) = -0.40247669763528    ;
%     n(27) =  0.58083399985759    ;
%     n(28) =  0.49969146990806E-2 ;
%     n(29) = -0.31358700712549E-1 ;
%     n(30) = -0.74315929710341    ;
%     n(31) =  0.47807329915480    ;
%     n(32) =  0.20527940895948E-1 ;
%     n(33) = -0.13636435110343    ;
%     n(34) =  0.14180634400617E-1 ;
%     n(35) =  0.83326504880713E-2 ;
%     n(36) = -0.29052336009585E-1 ;
%     n(37) =  0.38615085574206E-1 ;
%     n(38) = -0.20393486513704E-1 ;
%     n(39) = -0.16554050063734E-2 ;
%     n(40) =  0.19955571979541E-2 ;
%     n(41) =  0.15870308324157E-3 ;
%     n(42) = -0.16388568342530E-4 ;
%     n(43) =  0.43613615723811E-1 ;
%     n(44) =  0.34994005463765E-1 ;
%     n(45) = -0.76788197844621E-1 ;
%     n(46) =  0.22446277332006E-1 ;
%     n(47) = -0.62689710414685E-4 ;
%     n(48) = -0.55711118565645E-9 ;
%     n(49) = -0.19905718354408    ;
%     n(50) =  0.31777497330738    ;
%     n(51) = -0.11841182425981    ;
%     n(52) = -0.31306260323435E2  ;
%     n(53) =  0.31546140237781E2  ;
%     n(54) = -0.25213154341695E4  ;
%     n(55) = -0.14874640856724    ;
%     n(56) =  0.31806110878444    ;
% 
% end
%  
% function [alpha,beta,gamma,epsilon] = Coefficients_Group2()
%     
%     alpha(1)    = 20	;
%     alpha(2)    = 20	; 
%     alpha(3)    = 20	;
%     
%     beta(1)     = 150   ;
%     beta(2)     = 150   ;
%     beta(3)     = 250   ;
%     beta(4)     = 0.3   ;
%     beta(5)     = 0.3   ;
%     
%     gamma(1)    = 1.21	;
%     gamma(2)    = 1.21	;
%     gamma(3)    = 1.25	;
%     
%     epsilon(1)	= 1     ;
%     epsilon(2)	= 1     ;
%     epsilon(3)	= 1     ;
%     
% end
% 
% function [A,B,C,D,a,b] = Coefficients_Group3()
%    
%     A(1) = 0.32;
%     A(2) = 0.32;
%     
%     B(1) = 0.2;
%     B(2) = 0.2;
%     
%     C(1) = 28;
%     C(2) = 32;
%     
%     D(1) = 700;
%     D(2) = 800;
%     
%     a(1) = 3.5;
%     a(2) = 3.5;
%     
%     b(1) = 0.85;
%     b(2) = 0.95;
%     
% end
% 
