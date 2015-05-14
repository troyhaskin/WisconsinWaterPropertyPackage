function HelmDeriv = HelmholtzResidual_t(delta,tau)
    
    [c,d,t,n,alpha,beta,gamma,epsilon,A,B,C,D,a,b] = Coefficients_HelmholtzResidual();
    
    betaInv   = 1./beta         ;
	deltaMod  = (delta - 1).^2 + eps(delta)	;
    HelmDeriv = 0*delta               ;
    SumErr    = 0*delta               ; % Summation correction (Kahan summation)

    
    for k = 1:7
        Part      = n(k) * t(k) .* delta.^(d(k)) .* tau.^(t(k)-1)   ;
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
    
    for k = 8:51
        Part      = exp(-delta.^(c(k)))                                 ;
        Part      = n(k)*t(k) * delta.^(d(k)) .* tau.^(t(k)-1) .* Part  ;
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
    
    for k = 52:54
        m       = k - 51                                                        ;
        Arg     = -alpha(m)*(delta-epsilon(m)).^2 - beta(m)*(tau-gamma(m)).^2   ;
        Part    = t(k)./tau - 2.*beta(m).*(tau-gamma(m))                        ;
        Part    = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(Arg) .* Part       ;
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
    
    for k = 55:56
        m        = k - 51                                                        ;
        p        = k - 54                                                        ;

        Theta       = GetTheta    (deltaMod,tau,A(p),betaInv(m))   ;
        Delta       = GetDelta    (deltaMod,Theta,B(p),a(p))    ;
        Psi         = GetPsi      (deltaMod,tau,C(p),D(p))      ;
        Psi_t       = GetPsi_t    (tau,Psi,D(p))                ;
        Deltabi_t   = GetDeltabi_t(Delta,Theta,b(p))            ;
        
        Part      = n(k) * delta .*(Deltabi_t .* Psi + Delta.^(b(p)) .* Psi_t);
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
    
end

