function delG = EstimateDelGFromTau(tau)
    
    c    = [-2.03150240,-2.68302940,-5.38626492,-17.2991605,-44.7586581,-63.9201063];
    
    VarTheta = 1 - 1./tau;
    
    Part =        c(1) * VarTheta.^( 2/6);
    Part = Part + c(2) * VarTheta.^( 4/6);
    Part = Part + c(3) * VarTheta.^( 8/6);
    Part = Part + c(4) * VarTheta.^(18/6);
    Part = Part + c(5) * VarTheta.^(37/6);
    Part = Part + c(6) * VarTheta.^(71/6);
    
    delG = exp(Part);
    
end