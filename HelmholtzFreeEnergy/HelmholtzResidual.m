function Helm = HelmholtzResidual(delta,tau)
    
    [c,d,t,n,alpha,beta,gamma,epsilon,A,B,C,D,a,b] = Coefficients_HelmholtzResidual();
    
    deltaMod	= (delta - 1).^2 + eps(delta)    ;
    betaInv     = 1 ./ beta         ;
    Helm      = 0*delta               ;
    SumErr    = 0*delta               ; % Summation correction (Kahan summation)
    
    for k = 1:7
        Part    = n(k) * delta.^(d(k)) .* tau.^(t(k))   ;
        [Helm,SumErr] = KahanSum(Helm,Part,SumErr);
    end
    
    for k = 8:51
        Part    = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(-delta.^(c(k)))    ;
        [Helm,SumErr] = KahanSum(Helm,Part,SumErr);
    end
    
    for k = 52:54
        m       = k - 51                                                        ;
        Arg     = -alpha(m)*(delta-epsilon(m)).^2 - beta(m)*(tau-gamma(m)).^2   ;
        
        Part    = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(Arg)               ;
        [Helm,SumErr] = KahanSum(Helm,Part,SumErr);
    end
    
    for k = 55:56
        m        = k - 51	;
        p        = k - 54	;

        Theta       = GetTheta     (delta,tau,A(p),betaInv(m));
        Delta       = GetDelta     (deltaMod,Theta,B(p),a(p));
        Psi         = GetPsi       (deltaMod,tau,C(p),D(p));
        
        Part    = n(k) * Delta.^(b(p)) .* delta .* Psi  ;
        [Helm,SumErr] = KahanSum(Helm,Part,SumErr);
    end
    
end
