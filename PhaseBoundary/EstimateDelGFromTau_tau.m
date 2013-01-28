function delGdtau = EstimateDelGFromTau_tau(tau)
    
    c    = [-2.03150240,-2.68302940,-5.38626492,-17.2991605,-44.7586581,-63.9201063];
    
    VarTheta = 1 - 1./tau;
    
    Part =        c(1) * ( 2/6) * VarTheta.^(-4/6);
    Part = Part + c(2) * ( 4/6) * VarTheta.^(-2/6);
    Part = Part + c(3) * ( 8/6) * VarTheta.^( 2/6);
    Part = Part + c(4) * (18/6) * VarTheta.^(12/6);
    Part = Part + c(5) * (37/6) * VarTheta.^(31/6);
    Part = Part + c(6) * (71/6) * VarTheta.^(65/6);
    
    delG     = EstimateDelGFromTau(tau);
    delGdtau = delG .* Part ./ tau.^2;
end