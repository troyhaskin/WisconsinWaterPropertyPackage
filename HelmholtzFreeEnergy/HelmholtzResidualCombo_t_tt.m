function [Helm_t,Helm_tt] = HelmholtzResidualCombo_t_tt(delta,tau)
    
    persistent c d t n alpha beta gamma epsilon A B C D a b
    if isempty(c)
        [c,d,t,n,alpha,beta,gamma,epsilon,A,B,C,D,a,b] = ...
            Coefficients_HelmholtzResidual();
    end
    
    betaInv   = 1./beta         ;
	deltaMod  = (delta - 1).^2 + eps(delta)	;
    SumErr_t  = 0*delta               ; % Summation correction (Kahan summation)
    SumErr_tt = 0*delta               ; % Summation correction (Kahan summation)
    
%     for k = 1:7
%         Part                = n(k) * t(k) .* delta.^(d(k)) .* tau.^(t(k)-1)         ;
%         [Helm_t,SumErr_t]   = KahanSum(Helm_t,Part,SumErr_t)                        ;
%         [Helm_tt,SumErr_tt] = KahanSum(Helm_tt , Part./tau * (t(k)-1) , SumErr_tt)  ;
%     end
%
%   bsxfun version
     Part   = bsxfun(@times,n(1:7).*t(1:7),...
                bsxfun(@power,delta,d(1:7)).*bsxfun(@power,tau,t(1:7)-1));
    Helm_tt = sum(bsxfun(@times,bsxfun(@rdivide,Part,tau),t(1:7)-1),2);
    Helm_t  = sum(Part,2);
    

%     for k = 8:51
%         Part      = exp(-delta.^(c(k)))                                 ;
%         Part      = n(k)*t(k) * delta.^(d(k)) .* tau.^(t(k)-1) .* Part  ;
%         [Helm_t,SumErr_t] = KahanSum(Helm_t,Part,SumErr_t);
%         [Helm_tt,SumErr_tt] = KahanSum(Helm_tt , Part./tau * (t(k)-1) , SumErr_tt);
%     end
%
%   bsxfun version
    dd      = d(8:51);
    cc      = c(8:51);
    tt      = t(8:51);
    nn      = n(8:51);
    Part    = exp(-bsxfun(@power,delta,cc));
    Part    = bsxfun(@times,nn.*tt,bsxfun(@power,delta,dd) .* bsxfun(@power,tau,tt-1)) .* Part  ;
    Helm_t  = Helm_t  + sum(Part,2);
    Helm_tt = Helm_tt + sum(bsxfun(@times,bsxfun(@rdivide,Part,tau),tt-1),2);
    
    
    for k = 52:54
        m       = k - 51                                                        ;
        Arg     = -alpha(m)*(delta-epsilon(m)).^2 - beta(m)*(tau-gamma(m)).^2   ;
        Mult1   = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(Arg)               ;
        Mult2   = (t(k)./tau - 2.*beta(m).*(tau-gamma(m)))                      ;
        Part    = Mult1 .* Mult2                                                ;
        [Helm_t,SumErr_t] = KahanSum(Helm_t,Part,SumErr_t)                      ;


        Part    = Mult1 .* (Mult2.^2 - t(k)./tau.^2 - 2*beta(m))       ;
        [Helm_tt,SumErr_tt] = KahanSum(Helm_tt,Part,SumErr_tt)      ;
        
        
    end
    
    for k = 55:56
        m = k - 51  ;
        p = k - 54  ;

        Theta      = GetTheta     (deltaMod,tau,A(p),betaInv(m));
        Delta      = GetDelta     (deltaMod,Theta,B(p),a(p))    ;
        Psi        = GetPsi       (deltaMod,tau,C(p),D(p))      ;
        Psi_t      = GetPsi_t     (tau,Psi,D(p))                ;
        Psi_tt     = GetPsi_tt    (tau,Psi,D(p))                ;
        Deltabi_t  = GetDeltabi_t (Delta,Theta,b(p))            ;
        Deltabi_tt = GetDeltabi_tt(Delta,Theta,b(p))            ;
        
        Part              = n(k) * delta .*(Deltabi_t .* Psi + Delta.^(b(p)) .* Psi_t)  ;
        [Helm_t,SumErr_t] = KahanSum(Helm_t,Part,SumErr_t)                              ;


        Part                = Deltabi_tt .* Psi + 2* Deltabi_t .* Psi_t         ;
        Part                =  n(k) * delta .*(Part + Delta.^(b(p)) .* Psi_tt)  ;
        [Helm_tt,SumErr_tt] = KahanSum(Helm_tt,Part,SumErr_tt)                  ;
        
    end
    
end

